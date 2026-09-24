"""Offline tests only: no ENA or FASTQ requests in this test suite."""

import csv
import hashlib
import importlib.util
import io
from pathlib import Path
from urllib.error import URLError

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("fetch_run_manifest", ROOT / "scripts/fetch_run_manifest.py")
manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(manifest)


def record(files=2, **updates):
    result = dict.fromkeys(manifest.FIELDS, "")
    result.update(study_accession="PRJNA451237", secondary_study_accession="SRP142554",
                  run_accession="SRR7062125", experiment_accession="SRX3993086",
                  sample_accession="SAMN08965226", library_name="mRNA_F6-O",
                  experiment_alias="mRNA_F6-O",
                  library_layout="PAIRED", library_strategy="RNA-Seq",
                  fastq_ftp=";".join(f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR706/005/SRR7062125/read_{n}.fastq.gz" for n in range(files)),
                  fastq_md5=";".join(f"{n:032x}" for n in range(files)),
                  fastq_bytes=";".join(str(100 + n) for n in range(files)))
    result.update(updates)
    return result


def report(*rows, fields=manifest.FIELDS):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


@pytest.mark.parametrize("count", [0, 1, 2, 3, 5])
def test_all_file_counts_are_preserved_without_inferring_layout(count):
    payload = report(record(count))
    output = manifest.parse_report(payload, retrieved_utc="2026-09-24T00:00:00+00:00")[0]
    assert output["fastq_file_count"] == str(count)
    assert output["fastq_total_bytes"] == str(sum(100 + n for n in range(count)))
    assert output["fastq_ftp"] == record(count)["fastq_ftp"]
    assert output["fastq_md5"] == record(count)["fastq_md5"]
    assert output["library_layout"] == "PAIRED"
    assert output["source_report_sha256"] == hashlib.sha256(payload).hexdigest()
    assert output["fastq_https"].count("https://ftp.sra.ebi.ac.uk/") == count


@pytest.mark.parametrize("update", [
    {"fastq_md5": ""}, {"fastq_bytes": "100"}, {"fastq_md5": "a" * 31 + ";" + "b" * 32},
    {"fastq_md5": "g" * 32 + ";" + "b" * 32}, {"fastq_bytes": "100;0"},
    {"fastq_bytes": "100;-1"}, {"fastq_bytes": "100;1.0"}, {"fastq_bytes": "100;1e3"},
    {"fastq_bytes": "100; 3"}, {"fastq_bytes": "100;101;"}, {"fastq_bytes": "100;;101"},
    {"run_accession": "../../bad"}, {"experiment_accession": "not-an-accession"},
    {"sample_accession": "not-an-accession"}, {"study_accession": "PRJNA123", "secondary_study_accession": "SRP123"},
    {"library_layout": "paired"}, {"read_count": "1.5"}, {"base_count": "-1"},
])
def test_invalid_metadata_is_rejected(update):
    with pytest.raises(manifest.ManifestError):
        manifest.parse_report(report(record(**update)))


@pytest.mark.parametrize("location", [
    "http://ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk.evil.example/vol1/fastq/read.fastq.gz",
    "https://evil.ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz",
    "https://user@ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk:443/vol1/fastq/read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/../read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/%2e%2e/read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq//read.fastq.gz",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz?download=1",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz#x",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz\\bad",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/read\n.fastq.gz",
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/read.fastq.gz ",
    "file:///vol1/fastq/read.fastq.gz", "", "https://[broken",
])
def test_untrusted_fastq_urls_rejected(location):
    with pytest.raises(ValueError):
        manifest.fastq_https(location)


@pytest.mark.parametrize("prefix", ["", "ftp://", "https://"])
def test_trusted_fastq_url_normalization(prefix):
    suffix = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR706/005/SRR7062125/SRR7062125_1.fastq.gz"
    assert manifest.fastq_https(prefix + suffix) == "https://" + suffix


def test_duplicate_urls_rejected_after_normalization():
    first = record()["fastq_ftp"].split(";")[0]
    with pytest.raises(manifest.ManifestError, match="Duplicate FASTQ"):
        manifest.parse_report(report(record(fastq_ftp=first + ";https://" + first)))


@pytest.mark.parametrize("payload", [b"", b"<html>bad gateway</html>", b"run_accession\n", b"\xff"])
def test_missing_bad_or_wrong_encoding_schema_rejected(payload):
    with pytest.raises(manifest.ManifestError):
        manifest.parse_report(payload)


def test_empty_result_duplicate_header_duplicate_run_and_truncated_row_rejected():
    with pytest.raises(manifest.ManifestError, match="no runs"):
        manifest.parse_report(report())
    with pytest.raises(manifest.ManifestError, match="duplicate TSV"):
        manifest.parse_report(report(record(), fields=(*manifest.FIELDS, "run_accession")))
    with pytest.raises(manifest.ManifestError, match="Duplicate run"):
        manifest.parse_report(report(record(), record()))
    with pytest.raises(manifest.ManifestError, match="Malformed TSV row"):
        manifest.parse_report(report(record()).rsplit(b"\t", 1)[0] + b"\n")
    with pytest.raises(manifest.ManifestError, match="Malformed TSV row"):
        manifest.parse_report(report(record()).rstrip(b"\n") + b"\textra\n")


def test_reordered_headers_additional_columns_and_utf8_bom_are_supported():
    row = record(sample_title="Ovarian sample – metadata", future_ena_field="value")
    payload = b"\xef\xbb\xbf" + report(row, fields=(*reversed(manifest.FIELDS), "future_ena_field"))
    assert manifest.parse_report(payload)[0]["sample_title"] == row["sample_title"]


def test_controls_in_quoted_cells_rejected():
    with pytest.raises(manifest.ManifestError, match="Control character"):
        manifest.parse_report(report(record(sample_title="sample\nnewline")))


def test_expected_count_and_study_checked():
    with pytest.raises(manifest.ManifestError, match="Expected 31"):
        manifest.parse_report(report(record()), expected_runs=31)
    with pytest.raises(manifest.ManifestError, match="expected_runs"):
        manifest.parse_report(report(record()), expected_runs=0)
    assert manifest.parse_report(report(record()), study="PRJNA451237")[0]["run_accession"] == "SRR7062125"
    with pytest.raises(manifest.ManifestError):
        manifest.report_url("SRP142554&other=bad")


def test_round_trip_preserves_lists_and_sorting():
    rows = manifest.parse_report(report(record(run_accession="SRR7062126"), record()))
    actual = list(csv.DictReader(io.StringIO(manifest.manifest_tsv(rows)), delimiter="\t"))
    assert actual == rows
    assert actual[0]["run_accession"] == "SRR7062125"


def test_mapping_requires_alias_strategy_layout_and_paper_label_agreement():
    valid = manifest.parse_report(report(record()))[0]
    assert (valid["sample_id"], valid["assay"]) == ("F6-O", "mRNA")
    for update in [{"experiment_alias": "mRNA_F1-O"}, {"library_strategy": "miRNA-Seq"},
                   {"library_name": "mRNA_F9-O", "experiment_alias": "mRNA_F9-O"}]:
        row = manifest.parse_report(report(record(**update)))[0]
        assert row["sample_id"] == row["assay"] == ""
        assert row["mapping_evidence"] == "unresolved"
    for sample, expected in [("M2-EB", "M2-EB"), ("M2-EA", "")]:
        row = manifest.parse_report(report(record(1, library_name=sample, experiment_alias=sample,
                                                   library_layout="SINGLE", library_strategy="miRNA-Seq")))[0]
        assert row["sample_id"] == expected


def test_atomic_output_and_explicit_overwrite(tmp_path):
    output = tmp_path / "runs.tsv"
    manifest.atomic_write(output, "old\n")
    with pytest.raises(manifest.ManifestError, match="already exists"):
        manifest.atomic_write(output, "unapproved\n")
    assert output.read_text() == "old\n"
    manifest.atomic_write(output, "new\n", overwrite=True)
    assert output.read_text() == "new\n"
    assert list(tmp_path.iterdir()) == [output]


def test_failed_atomic_replacement_preserves_existing_file_and_cleans_temp(tmp_path, monkeypatch):
    output = tmp_path / "runs.tsv"
    output.write_text("old\n")
    def fail(*args):
        raise OSError("mock replace error")
    monkeypatch.setattr(manifest.os, "replace", fail)
    with pytest.raises(OSError):
        manifest.atomic_write(output, "new\n", overwrite=True)
    assert output.read_text() == "old\n"
    assert list(tmp_path.iterdir()) == [output]


def test_no_clobber_if_another_process_creates_destination(tmp_path, monkeypatch):
    output = tmp_path / "runs.tsv"
    link = manifest.os.link
    def create_first(source, target):
        output.write_text("concurrent\n")
        return link(source, target)
    monkeypatch.setattr(manifest.os, "link", create_first)
    with pytest.raises(FileExistsError):
        manifest.atomic_write(output, "new\n")
    assert output.read_text() == "concurrent\n"
    assert list(tmp_path.iterdir()) == [output]


def test_symlink_output_is_not_followed(tmp_path):
    original = tmp_path / "original.tsv"
    original.write_text("old\n")
    output = tmp_path / "runs.tsv"
    output.symlink_to(original)
    with pytest.raises(manifest.ManifestError, match="symlink"):
        manifest.atomic_write(output, "new\n", overwrite=True)
    assert original.read_text() == "old\n"


def test_offline_cli_never_uses_network_and_failure_preserves_output(tmp_path, monkeypatch):
    source = tmp_path / "report.tsv"
    source.write_bytes(report(record()))
    output = tmp_path / "runs.tsv"
    def fail(*args):
        pytest.fail("offline mode must not use network")
    monkeypatch.setattr(manifest, "fetch_report", fail)
    assert manifest.main(["--input-report", str(source), "--output", str(output)]) == 0
    previous = output.read_bytes()
    source.write_text("invalid\n")
    assert manifest.main(["--input-report", str(source), "--output", str(output), "--overwrite"]) == 1
    assert output.read_bytes() == previous


def test_fetch_failure_produces_no_output(tmp_path, monkeypatch):
    def fail(*args):
        raise URLError("offline")
    monkeypatch.setattr(manifest, "fetch_report", fail)
    output = tmp_path / "runs.tsv"
    assert manifest.main(["--output", str(output)]) == 1
    assert not output.exists()


def test_fetch_reads_metadata_only_with_size_limit(monkeypatch):
    seen = []
    class Response:
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def read(self, size):
            assert size == manifest.MAX_REPORT_BYTES + 1
            return report(record())
    class Opener:
        def open(self, request, timeout):
            seen.append(request.full_url)
            assert timeout == 30.0
            return Response()
    monkeypatch.setattr(manifest, "build_opener", lambda *args: Opener())
    manifest.fetch_report("SRP142554")
    assert len(seen) == 1 and seen[0].startswith(manifest.API + "?")
    assert "result=read_run" in seen[0]
    assert "fastq_bytes" in seen[0]


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf")])
def test_bad_timeout_rejected_before_network(timeout):
    with pytest.raises(manifest.ManifestError, match="timeout"):
        manifest.fetch_report("SRP142554", timeout)


def test_response_limit_and_untrusted_redirects(monkeypatch):
    monkeypatch.setattr(manifest, "MAX_REPORT_BYTES", 4)
    with pytest.raises(manifest.ManifestError, match="safety limit"):
        manifest.parse_report(b"12345")
    redirects = manifest.TrustedRedirects()
    for url in ["http://www.ebi.ac.uk/ena/portal/api/filereport", "https://evil.example/report", "https://www.ebi.ac.uk/other"]:
        with pytest.raises(manifest.ManifestError, match="redirected"):
            redirects.redirect_request(None, None, 302, "", {}, url)


def test_checked_in_snapshot_matches_paper_sample_assay_inventory():
    with (ROOT / "config/runs.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    validated = manifest.parse_report(
        report(*[{field: row[field] for field in manifest.FIELDS} for row in rows]),
        expected_runs=31,
    )
    expected = ({(sample, "mRNA") for sample in manifest.PAPER_SAMPLES}
                | {(sample, "miRNA") for sample in manifest.PAPER_SAMPLES - {"M2-EA"}})
    assert {(row["sample_id"], row["assay"]) for row in rows} == expected
    assert sum(int(row["fastq_file_count"]) for row in rows) == 47
    assert sum(int(row["fastq_total_bytes"]) for row in rows) == 145909488473
    assert len({row["source_report_sha256"] for row in rows}) == 1
    for actual, checked in zip(sorted(rows, key=lambda row: row["run_accession"]), validated):
        for field in ["sample_id", "assay", "mapping_evidence", "fastq_https", "fastq_file_count", "fastq_total_bytes"]:
            assert actual[field] == checked[field]
