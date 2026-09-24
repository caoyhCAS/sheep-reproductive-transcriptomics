#!/usr/bin/env python3
"""Plan or explicitly execute a method-based reconstruction of Heredity 2018.

The default action is read-only JSON planning. No downloads or tool installation.
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


DOI = "10.1038/s41437-018-0090-1"
STAGES = ("mrna", "cuffdiff", "trinity", "rsem", "denovo", "mirna", "all")
EXPECTED = {}
for individual in range(1, 7):
    EXPECTED[f"F{individual}-O"] = ("Finnsheep", "ovary", f"F{individual}", "single")
for individual in (2, 3):
    EXPECTED[f"F{individual}-E"] = ("Finnsheep", "endometrium", f"F{individual}", "single")
for individual in (1, 2, 3):
    for replicate in ("A", "B"):
        EXPECTED[f"M{individual}-O{replicate}"] = ("European_mouflon", "ovary", f"M{individual}", replicate)
for replicate in ("A", "B"):
    EXPECTED[f"M2-E{replicate}"] = ("European_mouflon", "endometrium", "M2", replicate)
LEGACY = {"fastqc": "0.11.4", "cutadapt": "1.9.1", "tophat": "2.0.8b",
          "cufflinks": "2.1.1", "cuffmerge": "2.1.1", "cuffdiff": "2.1.1",
          "trinity": "2.1.1", "prinseq": "0.20.4"}


def valid_path(value, label):
    # Some archived tools internally build shell strings. Safe outer argv alone
    # is insufficient, so paths intentionally exclude whitespace/metacharacters.
    if (not isinstance(value, str) or not value or value.startswith("-") or
            not (re.fullmatch(r"[A-Za-z0-9_./-]+", value) or re.fullmatch(r"<UNRESOLVED:[A-Za-z0-9_.]+>", value))):
        raise ValueError(f"{label}: use only ASCII letters, digits, underscore, dot, slash or hyphen in legacy paths")
    return value


def read_samples(path):
    columns = {"sample_id", "species", "tissue", "individual_id", "tissue_replicate", "mrna_r1", "mrna_r2", "mirna_fastq"}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not columns <= set(reader.fieldnames or []):
            raise ValueError(f"Sample table requires columns: {sorted(columns)}")
        rows = list(reader)
    if len(rows) != 16 or {r["sample_id"] for r in rows} != set(EXPECTED):
        raise ValueError("Paper manifest must contain exactly the 16 unique reported mRNA sample IDs")
    paths = []
    for row in rows:
        sample = row["sample_id"]
        actual = tuple(row[c] for c in ("species", "tissue", "individual_id", "tissue_replicate"))
        if actual != EXPECTED[sample]:
            raise ValueError(f"Incorrect biological individual/tissue/replicate metadata for {sample}")
        for key in ("mrna_r1", "mrna_r2"):
            paths.append(valid_path(row[key], f"{sample}.{key}"))
        if sample == "M2-EA":
            if row["mirna_fastq"] != ".":
                raise ValueError("M2-EA has no miRNA library; its mirna_fastq must be '.'")
        else:
            paths.append(valid_path(row["mirna_fastq"], f"{sample}.mirna_fastq"))
    if len(paths) != len({Path(path).resolve() for path in paths}):
        raise ValueError("Input FASTQ paths must be unique across mates/libraries")
    return rows


def plan(config, samples, stage="all"):
    if stage not in STAGES:
        raise ValueError(f"Unknown stage {stage}")
    root = valid_path(config.get("output_dir"), "output_dir")
    if Path(root).resolve() == Path.cwd().resolve() or Path(root).resolve() == Path("/"):
        raise ValueError("Use a dedicated output directory")
    threads = config.get("threads")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    steps, blockers, acknowledgments = [], [], {"method_reconstruction", "legacy_environment", "defaults_unreported"}
    tools_used = set()
    paths = lambda *parts: str(Path(root).joinpath(*parts))
    selected = {stage} if stage != "all" else {"mrna", "cuffdiff", "denovo", "mirna"}
    do_diff = "cuffdiff" in selected
    do_map = "mrna" in selected or do_diff
    do_trinity = bool({"trinity", "denovo"} & selected)
    do_rsem = bool({"rsem", "denovo"} & selected)
    do_mirna = "mirna" in selected

    def require(section, key, choices=None):
        value = config.get(section, {}).get(key)
        if value is None or value == "" or (choices is not None and value not in choices):
            blockers.append(f"Set {section}.{key}" + (f" to one of {choices}" if choices else ""))
            return f"<UNRESOLVED:{section}.{key}>"
        return value

    def command(tool):
        tools_used.add(tool)
        value = config.get("tools", {}).get(tool, {}).get("command", [tool])
        if not isinstance(value, list) or not value or not all(isinstance(v, str) and v and "\x00" not in v for v in value):
            raise ValueError(f"tools.{tool}.command must be a nonempty argument list")
        return value

    def add(identifier, stage_name, tool, args, inputs, outputs, depends=(), files=None, stdout=None, metadata=None):
        for path in [*inputs, *outputs, *list((files or {}).keys())]:
            valid_path(path, identifier)
        step = {"id": identifier, "stage": stage_name, "argv": command(tool) + [str(v) for v in args],
                "inputs": inputs, "outputs": outputs, "depends_on": list(depends),
                "generated_files": files or {}, "stdout": stdout, "metadata": metadata or {},
                "log": paths("logs", identifier + ".log")}
        steps.append(step)

    def adapter(section, key):
        value = require(section, key)
        if not value.startswith("<UNRESOLVED:") and not re.fullmatch(r"[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+", value):
            raise ValueError(f"{section}.{key}: provide the actual IUPAC adapter sequence")
        return value

    def fq_stem(path):
        name = Path(path).name
        for suffix in (".gz", ".fastq", ".fq"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        return name

    mrna_rows = samples if do_map else ([r for r in samples if r["species"] == "European_mouflon"] if do_trinity or do_rsem else [])
    trimmed = {}
    if mrna_rows:
        a1, a2 = adapter("mrna", "adapter_r1"), adapter("mrna", "adapter_r2")
        require("mrna", "adapter_resolution_note")
        for row in mrna_rows:
            sample = row["sample_id"]
            qc = paths("mrna", "fastqc", sample)
            raw = [row["mrna_r1"], row["mrna_r2"]]
            qc_out = [str(Path(qc) / (fq_stem(p) + "_fastqc.html")) for p in raw]
            if len(set(qc_out)) != 2:
                raise ValueError(f"{sample}: mate filenames must have distinct basenames for FastQC")
            add(f"mrna_fastqc_{sample}", "mrna", "fastqc", ["--outdir", qc, *raw], raw, qc_out)
            reads = [paths("mrna", "trimmed", sample + ".R1.fastq"), paths("mrna", "trimmed", sample + ".R2.fastq")]
            add(f"mrna_trim_{sample}", "mrna", "cutadapt", ["-a", a1, "-A", a2, "-o", reads[0], "-p", reads[1], *raw], raw, reads, [f"mrna_fastqc_{sample}"])
            trimmed[sample] = reads

    groups = {}
    for row in samples:
        groups.setdefault(row["species"] + "_" + row["tissue"], []).append(row["sample_id"])
    bams, assemblies = {}, {}
    if do_map:
        acknowledgments.add("reference_annotation_conflict")
        genome, gtf = require("reference", "genome_fasta"), require("reference", "annotation_gtf")
        assembly = require("reference", "genome_assembly")
        gtf_assembly = require("reference", "annotation_assembly")
        require("reference", "resolution_note")
        if not gtf_assembly.startswith("<UNRESOLVED:") and assembly != gtf_assembly:
            blockers.append("Genome and GTF assemblies differ; provide a consistent pair, not just acknowledgment")
        prefix = require("reference", "bowtie2_index_prefix")
        index_files = config.get("reference", {}).get("bowtie2_index_files")
        if not isinstance(index_files, list) or len(index_files) != 6:
            blockers.append("reference.bowtie2_index_files must explicitly name all six matching index files")
            index_files = [prefix + suffix for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")]
        elif set(index_files) not in [{prefix + suffix + extension for suffix in (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")} for extension in (".bt2", ".bt2l")]:
            blockers.append("reference.bowtie2_index_files do not match the six files for bowtie2_index_prefix")
        library = require("mrna", "library_type", ("fr-unstranded", "fr-firststrand", "fr-secondstrand"))
        guide = require("mrna", "cufflinks_annotation_mode", ("guided", "unguided"))
        for row in samples:
            sample = row["sample_id"]
            out = paths("mrna", "tophat", sample)
            bam = str(Path(out) / "accepted_hits.bam")
            add(f"tophat_{sample}", "mrna", "tophat", ["-p", threads, "--library-type", library, "-G", gtf, "-o", out, prefix, *trimmed[sample]], [*trimmed[sample], genome, gtf, *index_files], [bam], [f"mrna_trim_{sample}"])
            out = paths("mrna", "cufflinks", sample)
            transcript = str(Path(out) / "transcripts.gtf")
            annotation_args = ["-g", gtf] if guide == "guided" else []
            add(f"cufflinks_{sample}", "mrna", "cufflinks", ["-p", threads, "--library-type", library, *annotation_args, "-o", out, bam], [bam, *([gtf] if guide == "guided" else [])], [transcript], [f"tophat_{sample}"])
            bams[sample], assemblies[sample] = bam, transcript

    contrasts = {
        "ovary_species": ("European_mouflon_ovary", "Finnsheep_ovary"),
        "endometrium_species": ("European_mouflon_endometrium", "Finnsheep_endometrium"),
        "mouflon_tissues": ("European_mouflon_ovary", "European_mouflon_endometrium"),
    }
    if do_diff:
        acknowledgments |= {"pseudoreplication", "merged_gtf_ambiguity"}
        policy = require("cuffdiff", "transcriptome_policy", ("union_per_contrast", "provided_per_contrast"))
        require("cuffdiff", "resolution_note")
        merged = {}
        for group, members in groups.items():
            directory = paths("mrna", "cuffmerge", group)
            listing = paths("mrna", "assembly_lists", group + ".txt")
            output = str(Path(directory) / "merged.gtf")
            add(f"cuffmerge_{group}", "cuffdiff", "cuffmerge", ["-p", threads, "-g", gtf, "-s", genome, "-o", directory, listing], [genome, gtf, *[assemblies[s] for s in members]], [output], [f"cufflinks_{s}" for s in members], {listing: "\n".join(assemblies[s] for s in members) + "\n"})
            merged[group] = output
        for name, (numerator, denominator) in contrasts.items():
            deps = [f"cuffmerge_{numerator}", f"cuffmerge_{denominator}"]
            if policy == "provided_per_contrast":
                transcript = config.get("cuffdiff", {}).get("contrast_gtfs", {}).get(name)
                if not transcript:
                    blockers.append(f"Set cuffdiff.contrast_gtfs.{name}")
                    transcript = f"<UNRESOLVED:cuffdiff.contrast_gtfs.{name}>"
            else:
                directory = paths("mrna", "contrast_gtf", name)
                transcript = str(Path(directory) / "merged.gtf")
                listing = paths("mrna", "assembly_lists", name + ".contrast.txt")
                add(f"merge_contrast_{name}", "cuffdiff", "cuffmerge", ["-p", threads, "-g", gtf, "-s", genome, "-o", directory, listing], [genome, gtf, merged[numerator], merged[denominator]], [transcript], deps, {listing: merged[denominator] + "\n" + merged[numerator] + "\n"})
                deps = [f"merge_contrast_{name}"]
            directory = paths("mrna", "cuffdiff", name)
            denominator_bams = [bams[s] for s in groups[denominator]]
            numerator_bams = [bams[s] for s in groups[numerator]]
            metadata = {"numerator": numerator, "denominator": denominator,
                        "sample_1": denominator, "sample_2": numerator,
                        "log2fc_definition": "log2(value_2/value_1)", "biological_replicates": {g: sorted({EXPECTED[s][2] for s in groups[g]}) for g in (numerator, denominator)}}
            add(f"cuffdiff_{name}", "cuffdiff", "cuffdiff", ["-p", threads, "--library-type", library, "-L", denominator + "," + numerator, "-o", directory, transcript, ",".join(denominator_bams), ",".join(numerator_bams)], [transcript, *denominator_bams, *numerator_bams], [str(Path(directory) / "gene_exp.diff")], deps, metadata=metadata)

    if do_trinity:
        strand = require("trinity", "strand", ("unstranded", "RF", "FR"))
        memory = config.get("trinity", {}).get("max_memory", "32G")
        if not isinstance(memory, str) or not re.fullmatch(r"[1-9][0-9]*G", memory):
            raise ValueError("trinity.max_memory must be e.g. '32G'")
        members = [r["sample_id"] for r in samples if r["species"] == "European_mouflon"]
        directory = paths("denovo", "trinity")
        trinity_fasta = str(Path(directory) / "Trinity.fasta")
        strand_args = [] if strand == "unstranded" else ["--SS_lib_type", strand]
        add("trinity_all_eight_mouflon", "denovo", "trinity", ["--seqType", "fq", "--no_version_check", "--max_memory", memory, "--CPU", threads, "--left", ",".join(trimmed[s][0] for s in members), "--right", ",".join(trimmed[s][1] for s in members), *strand_args, "--output", directory], [p for s in members for p in trimmed[s]], [trinity_fasta], [f"mrna_trim_{s}" for s in members], metadata={"sample_ids": members})
    else:
        trinity_fasta = config.get("rsem", {}).get("trinity_fasta", paths("denovo", "trinity", "Trinity.fasta"))
    if do_rsem:
        acknowledgments |= {"rsem_unreported_settings", "pseudoreplication"}
        aligner = require("rsem", "aligner", ("bowtie", "bowtie2"))
        strand = require("rsem", "strandedness", ("none", "forward", "reverse"))
        require("rsem", "resolution_note")
        map_path = paths("denovo", "rsem_reference", "gene_transcript_map.tsv")
        deps = ["trinity_all_eight_mouflon"] if do_trinity else []
        add("rsem_gene_map", "rsem", "rsem_extract", [trinity_fasta, map_path], [trinity_fasta], [map_path], deps)
        prefix = paths("denovo", "rsem_reference", "mouflon")
        add("rsem_reference", "rsem", "rsem_prepare", ["--transcript-to-gene-map", map_path, "--" + aligner, trinity_fasta, prefix], [map_path, trinity_fasta], [prefix + ".grp", prefix + ".ti", prefix + ".seq"], ["rsem_gene_map"])
        for sample in ("M2-OA", "M2-OB", "M2-EA", "M2-EB"):
            prefix_out = paths("denovo", "rsem", sample)
            # Bowtie is the RSEM default; only Bowtie2 is an explicit quantification flag.
            align_args = ["--bowtie2"] if aligner == "bowtie2" else []
            add("rsem_" + sample, "rsem", "rsem_quant", ["--paired-end", "-p", threads, "--strandedness", strand, *align_args, *trimmed[sample], prefix, prefix_out], [*trimmed[sample], prefix + ".grp", prefix + ".ti", prefix + ".seq"], [prefix_out + ".genes.results", prefix_out + ".isoforms.results"], ["rsem_reference", f"mrna_trim_{sample}"], metadata={"individual_id": "M2", "tissue": EXPECTED[sample][1]})

    if do_mirna:
        acknowledgments.add("mirna_alignment_conflict")
        adapter_mirna = adapter("mirna", "adapter")
        require("mirna", "adapter_resolution_note")
        for row in samples:
            if row["mirna_fastq"] == ".":
                continue
            sample, raw = row["sample_id"], row["mirna_fastq"]
            directory = paths("mirna", "fastqc", sample)
            add("mirna_fastqc_" + sample, "mirna", "fastqc", ["--outdir", directory, raw], [raw], [str(Path(directory) / (fq_stem(raw) + "_fastqc.html"))])
            trimmed_mirna = paths("mirna", "trimmed", sample + ".fastq")
            add("mirna_trim_" + sample, "mirna", "cutadapt", ["-a", adapter_mirna, "-o", trimmed_mirna, raw], [raw], [trimmed_mirna], ["mirna_fastqc_" + sample])
            prefix = paths("mirna", "clean", sample)
            add("prinseq_" + sample, "mirna", "prinseq", ["-fastq", trimmed_mirna, "-min_len", "18", "-max_len", "26", "-ns_max_n", "0", "-noniupac", "-out_format", "3", "-out_good", prefix, "-out_bad", "null"], [trimmed_mirna], [prefix + ".fastq"], ["mirna_trim_" + sample])

    for tool in sorted(tools_used):
        spec = config.get("tools", {}).get(tool, {})
        if tool in LEGACY and spec.get("expected_version") != LEGACY[tool]:
            blockers.append(f"tools.{tool}.expected_version must remain reported version {LEGACY[tool]}")
        if tool in ("rsem_quant", "rsem_prepare", "rsem_extract") and not spec.get("expected_version"):
            blockers.append(f"Set tools.{tool}.expected_version: RSEM version was not reported")
        if not spec.get("version_command"):
            blockers.append(f"Set tools.{tool}.version_command for installed-environment verification")
    for acknowledgment in sorted(acknowledgments):
        if config.get("acknowledgments", {}).get(acknowledgment) is not True:
            blockers.append(f"Acknowledge {acknowledgment} in config after reviewing docs/WORKFLOW.md")
    return {"doi": DOI, "status": "method_reconstruction_not_original_code", "stage": stage,
            "library_counts": {"mrna": 16, "mirna": 15, "biological_individuals": 9},
            "execution_blockers": list(dict.fromkeys(blockers)), "required_acknowledgments": sorted(acknowledgments),
            "steps": steps, "tools_used": sorted(tools_used),
            "handoffs": ["miRDeep2 mapper.pl/miRDeep2.pl, alignment format and Bowtie/Oar assembly conflict require a reviewed separate mapping plan; not executed here.",
                         "DESeq1.22.0, edgeR, BLAST, PANTHER, DAVID, TargetScan and Cytoscape are not executed by this planner.",
                         "Mouflon tissue samples are repeated tissues from three individuals; endometrium comes only from M2. RSEM four samples all come from M2."]}


def execution_preflight(plan_data, config):
    valid_path(str(Path.cwd()), "working directory")
    errors = list(plan_data["execution_blockers"])
    output_root = Path(config["output_dir"])
    if output_root.exists() and (not output_root.is_dir() or any(output_root.iterdir())):
        errors.append("output_dir already contains files; choose a new directory to protect existing tool outputs")
    produced = {p for step in plan_data["steps"] for p in step["outputs"]}
    generated = {p for step in plan_data["steps"] for p in step["generated_files"]}
    for path in sorted({p for step in plan_data["steps"] for p in step["inputs"]} - produced - generated):
        if not Path(path).is_file() or Path(path).stat().st_size == 0:
            errors.append(f"Missing or empty external input: {path}")
    for path in produced | generated:
        if Path(path).exists():
            errors.append(f"Refusing to overwrite output: {path}; choose a new output_dir")
    if errors:
        raise ValueError("Execution blocked:\n- " + "\n- ".join(errors))
    versions = {}
    for tool in plan_data["tools_used"]:
        spec = config["tools"][tool]
        for argv in (spec["command"], spec["version_command"]):
            if not isinstance(argv, list) or not argv or not all(isinstance(a, str) and a for a in argv):
                raise ValueError(f"Invalid command argument list for {tool}")
            if shutil.which(argv[0]) is None:
                raise ValueError(f"Executable unavailable: {argv[0]}")
        result = subprocess.run(spec["version_command"], text=True, capture_output=True, timeout=30, check=False)
        version = (result.stdout + result.stderr).strip()
        expected = spec["expected_version"]
        if result.returncode not in spec.get("version_exit_codes", [0]) or re.search(r"(?<![A-Za-z0-9.])v?" + re.escape(expected) + r"(?![A-Za-z0-9.])", version) is None:
            raise ValueError(f"Version check failed for {tool}: expected {expected}, received {version[:300]!r}")
        versions[tool] = {"expected": expected, "reported": version}
        if tool == "cuffmerge":
            self_command = spec.get("self_version_command")
            if not isinstance(self_command, list) or not self_command or not all(isinstance(a, str) and a for a in self_command):
                raise ValueError("Set tools.cuffmerge.self_version_command to verify the separate script version")
            result = subprocess.run(self_command, text=True, capture_output=True, timeout=30, check=False)
            own_version = (result.stdout + result.stderr).strip()
            own_expected = spec.get("self_expected_version", "1.0.0")
            if result.returncode != 0 or re.search(r"(?<![A-Za-z0-9.])v?" + re.escape(own_expected) + r"(?![A-Za-z0-9.])", own_version) is None:
                raise ValueError("Cuffmerge script version check failed")
            versions[tool]["script_version"] = own_version
    return versions


def execute(plan_data, config):
    versions = execution_preflight(plan_data, config)
    root = Path(config["output_dir"])
    root.mkdir(parents=True, exist_ok=True)
    (root / "execution_plan.json").write_text(json.dumps(plan_data, indent=2) + "\n", encoding="utf-8")
    (root / "verified_versions.json").write_text(json.dumps(versions, indent=2) + "\n", encoding="utf-8")
    completed = set()
    for step in plan_data["steps"]:
        if not set(step["depends_on"]) <= completed:
            raise ValueError(f"Unmet dependency in {step['id']}")
        for path in step["inputs"]:
            if not Path(path).is_file():
                raise ValueError(f"Input was not generated for {step['id']}: {path}")
        for path in [*step["outputs"], *step["generated_files"], step["log"]]:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
        for path, content in step["generated_files"].items():
            Path(path).write_text(content, encoding="utf-8")
        # FastQC requires its explicit output directory to exist.
        if "fastqc" in step["id"]:
            Path(step["outputs"][0]).parent.mkdir(parents=True, exist_ok=True)
        with open(step["log"], "w", encoding="utf-8") as log:
            subprocess.run(step["argv"], stdout=log, stderr=subprocess.STDOUT, check=True)
        missing = [p for p in step["outputs"] if not Path(p).is_file()]
        if missing:
            raise ValueError(f"{step['id']} finished without expected outputs: {missing}")
        completed.add(step["id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/paper.json")
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--stage", choices=STAGES, default="all")
    parser.add_argument("--output", help="Write a reviewable JSON plan instead of printing it")
    parser.add_argument("--execute", action="store_true", help="Run selected legacy commands only after all preflight checks pass")
    args = parser.parse_args()
    try:
        with open(args.config, encoding="utf-8") as handle:
            config = json.load(handle)
        result = plan(config, read_samples(args.samples), args.stage)
        if args.output:
            destination = Path(args.output)
            if destination.resolve() in {Path(args.config).resolve(), Path(args.samples).resolve()}:
                raise ValueError("Plan output must not replace the input config or sample manifest")
            if destination.exists():
                raise ValueError("Plan output already exists; choose a new path to preserve prior plans")
            if args.execute and (destination.resolve() == Path(config["output_dir"]).resolve() or Path(config["output_dir"]).resolve() in destination.resolve().parents):
                raise ValueError("With --execute, put --output outside output_dir; execution writes its own plan")
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        else:
            print(json.dumps(result, indent=2))
        if args.execute:
            execute(result, config)
    except (ValueError, TypeError, OSError, subprocess.SubprocessError) as exc:
        parser.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
