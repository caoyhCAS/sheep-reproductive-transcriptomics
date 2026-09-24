"""Planner/guard tests, not biological validation of archived executables."""
import copy
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("workflow", ROOT / "scripts" / "workflow.py")
WORKFLOW = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WORKFLOW)


def inputs():
    return json.loads((ROOT / "config" / "paper.json").read_text()), WORKFLOW.read_samples(ROOT / "config" / "samples.tsv")


def test_exact_library_and_biological_sample_counts():
    config, rows = inputs()
    assert len(rows) == 16
    assert sum(r["mirna_fastq"] != "." for r in rows) == 15
    assert len({r["individual_id"] for r in rows}) == 9
    assert {r["individual_id"] for r in rows if r["species"] == "European_mouflon" and r["tissue"] == "endometrium"} == {"M2"}


def write_manifest(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_manifest_rejects_alias_of_same_fastq(tmp_path):
    _, rows = inputs()
    rows[1]["mrna_r1"] = "./" + rows[0]["mrna_r1"]
    path = tmp_path / "samples.tsv"
    write_manifest(path, rows)
    with pytest.raises(ValueError, match="unique"):
        WORKFLOW.read_samples(path)


@pytest.mark.parametrize("change", ["invent_mirna", "wrong_individual", "duplicate_sample"])
def test_manifest_rejects_pseudoreplication_metadata_errors(tmp_path, change):
    _, rows = inputs()
    if change == "invent_mirna":
        next(r for r in rows if r["sample_id"] == "M2-EA")["mirna_fastq"] = "invented.fastq"
    elif change == "wrong_individual":
        next(r for r in rows if r["sample_id"] == "M2-EB")["individual_id"] = "M3"
    else:
        rows[1] = rows[0].copy()
    path = tmp_path / "samples.tsv"
    write_manifest(path, rows)
    with pytest.raises(ValueError):
        WORKFLOW.read_samples(path)


def test_default_plan_exposes_unresolved_choices_without_subprocess(tmp_path):
    config, rows = inputs()
    config["output_dir"] = str(tmp_path / "outputs")
    with patch.object(WORKFLOW.subprocess, "run", side_effect=AssertionError("default plan must not invoke tools")):
        result = WORKFLOW.plan(config, rows)
    assert not Path(config["output_dir"]).exists()
    assert any("reference.annotation_assembly" in x for x in result["execution_blockers"])
    assert any("mrna.adapter_r1" in x for x in result["execution_blockers"])
    assert any("cuffdiff.transcriptome_policy" in x for x in result["execution_blockers"])
    assert any("pseudoreplication" in x for x in result["execution_blockers"])


def test_full_plan_dependencies_trinity_and_rsem_members():
    config, rows = inputs()
    result = WORKFLOW.plan(config, rows)
    completed = set()
    for step in result["steps"]:
        assert set(step["depends_on"]) <= completed
        assert step["id"] not in completed
        assert isinstance(step["argv"], list)
        completed.add(step["id"])
    trinity = next(s for s in result["steps"] if s["id"] == "trinity_all_eight_mouflon")
    assert len(trinity["metadata"]["sample_ids"]) == 8
    assert len(trinity["argv"][trinity["argv"].index("--left") + 1].split(",")) == 8
    assert "--no_version_check" in trinity["argv"]
    rsem = [s for s in result["steps"] if s["id"] in {"rsem_M2-OA", "rsem_M2-OB", "rsem_M2-EA", "rsem_M2-EB"}]
    assert len(rsem) == 4
    assert all(s["metadata"]["individual_id"] == "M2" for s in rsem)
    mirna = [s for s in result["steps"] if s["id"].startswith("prinseq_")]
    assert len(mirna) == 15
    assert not any("M2-EA" in s["id"] for s in mirna)
    for step in mirna:
        args = step["argv"]
        assert args[args.index("-min_len") + 1] == "18"
        assert args[args.index("-max_len") + 1] == "26"
        assert args[args.index("-ns_max_n") + 1] == "0"
        assert "-noniupac" in args


def test_cuffdiff_direction_keeps_value2_numerator():
    config, rows = inputs()
    result = WORKFLOW.plan(config, rows, "cuffdiff")
    diffs = [s for s in result["steps"] if s["id"].startswith("cuffdiff_")]
    assert len(diffs) == 3
    for step in diffs:
        args, metadata = step["argv"], step["metadata"]
        assert args[args.index("-L") + 1] == metadata["denominator"] + "," + metadata["numerator"]
        assert metadata["sample_2"] == metadata["numerator"]
    endometrium = next(s for s in diffs if s["id"] == "cuffdiff_endometrium_species")
    assert endometrium["metadata"]["biological_replicates"]["European_mouflon_endometrium"] == ["M2"]


def test_preflight_blocks_before_tools_or_output_mutation(tmp_path):
    config, rows = inputs()
    config["output_dir"] = str(tmp_path / "output")
    result = WORKFLOW.plan(config, rows, "mrna")
    with patch.object(WORKFLOW.subprocess, "run", side_effect=AssertionError("blocked execution invoked a tool")):
        with pytest.raises(ValueError, match="Execution blocked"):
            WORKFLOW.execution_preflight(result, config)
    assert not Path(config["output_dir"]).exists()


@pytest.mark.parametrize("path", ["data/file name.fastq", "data/a;touch_PWN.fastq", "data/a$(id).fastq", "data/a`id`.fastq", "data/a,b.fastq"])
def test_reject_paths_unsafe_in_internal_legacy_shells(path):
    with pytest.raises(ValueError, match="ASCII"):
        WORKFLOW.valid_path(path, "input")


def test_reference_conflict_cannot_be_acknowledged_away():
    config, rows = inputs()
    config["reference"]["annotation_assembly"] = "Oar_v3.1"
    config["reference"]["resolution_note"] = "Intentional mismatch is still invalid"
    config["acknowledgments"]["reference_annotation_conflict"] = True
    result = WORKFLOW.plan(config, rows, "mrna")
    assert any("Genome and GTF assemblies differ" in x for x in result["execution_blockers"])


def test_version_probe_allows_trinity_v_prefix_and_legacy_exit_one(tmp_path):
    config = {"output_dir": str(tmp_path / "output"), "tools": {"trinity": {"command": ["Trinity"], "version_command": ["Trinity", "--version"], "expected_version": "2.1.1", "version_exit_codes": [0, 1]}}}
    plan_data = {"execution_blockers": [], "tools_used": ["trinity"], "steps": []}
    with patch.object(WORKFLOW.shutil, "which", return_value="/usr/bin/Trinity"), patch.object(WORKFLOW.subprocess, "run", return_value=SimpleNamespace(returncode=1, stdout="Trinity version: v2.1.1\n", stderr="")):
        report = WORKFLOW.execution_preflight(plan_data, config)
    assert report["trinity"]["expected"] == "2.1.1"


def test_cuffmerge_bundle_and_script_versions_are_separate(tmp_path):
    config, _ = inputs()
    config["output_dir"] = str(tmp_path / "outputs")
    plan_data = {"execution_blockers": [], "tools_used": ["cuffmerge"], "steps": []}
    responses = [SimpleNamespace(returncode=0, stdout="cufflinks v2.1.1", stderr=""), SimpleNamespace(returncode=0, stdout="merge_cuff_asms v1.0.0", stderr="")]
    with patch.object(WORKFLOW.shutil, "which", return_value="/usr/bin/tool"), patch.object(WORKFLOW.subprocess, "run", side_effect=responses):
        report = WORKFLOW.execution_preflight(plan_data, config)
    assert report["cuffmerge"]["script_version"] == "merge_cuff_asms v1.0.0"


def test_cli_cannot_overwrite_input_configuration(tmp_path):
    config = tmp_path / "config.json"
    original = (ROOT / "config" / "paper.json").read_text()
    config.write_text(original)
    result = subprocess.run([sys.executable, str(ROOT / "scripts" / "workflow.py"), "--config", str(config), "--samples", str(ROOT / "config" / "samples.tsv"), "--output", str(config)], capture_output=True, text=True)
    assert result.returncode == 2
    assert "must not replace" in result.stderr
    assert config.read_text() == original
