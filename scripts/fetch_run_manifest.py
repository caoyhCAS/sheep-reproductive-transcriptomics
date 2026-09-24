#!/usr/bin/env python3
"""Fetch ENA metadata only; never download FASTQ or infer biological conditions."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from urllib.parse import urlencode, urlsplit, urlunsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener


API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = (
    "study_accession", "secondary_study_accession", "run_accession",
    "experiment_accession", "sample_accession", "secondary_sample_accession",
    "experiment_alias", "library_name", "sample_alias", "sample_title",
    "library_layout", "library_strategy", "library_source", "library_selection",
    "instrument_model", "read_count", "base_count", "fastq_ftp", "fastq_md5",
    "fastq_bytes",
)
OUTPUT_FIELDS = (*FIELDS, "sample_id", "assay", "mapping_evidence",
                 "fastq_https", "fastq_file_count", "fastq_total_bytes",
                 "source_url", "retrieved_utc", "source_report_sha256")
MAX_REPORT_BYTES = 16 * 1024 * 1024
STUDY_PATTERN = r"(?:[SED]RP\d+|PRJ(?:NA|EB|DB)\d+)"
PAPER_SAMPLES = ({f"F{n}-O" for n in range(1, 7)} | {"F2-E", "F3-E"}
                 | {f"M{n}-O{rep}" for n in range(1, 4) for rep in "AB"}
                 | {"M2-EA", "M2-EB"})


class ManifestError(ValueError):
    """The report cannot safely be used as a metadata manifest."""


def report_url(study):
    if not isinstance(study, str) or not re.fullmatch(STUDY_PATTERN, study):
        raise ManifestError("study must be an INSDC study/BioProject accession")
    return API + "?" + urlencode({
        "accession": study, "result": "read_run", "fields": ",".join(FIELDS),
        "format": "tsv", "download": "false",
    })


class TrustedRedirects(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        target = urlsplit(newurl)
        if (target.scheme != "https" or target.netloc != "www.ebi.ac.uk"
                or target.path != "/ena/portal/api/filereport"):
            raise ManifestError("ENA report redirected outside the trusted HTTPS endpoint")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def fetch_report(study, timeout=30.0):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ManifestError("timeout must be finite and positive")
    request = Request(report_url(study), headers={
        "Accept": "text/plain", "User-Agent": "sheep-reproductive-transcriptomics/1.0",
    })
    with build_opener(TrustedRedirects()).open(request, timeout=timeout) as response:
        payload = response.read(MAX_REPORT_BYTES + 1)
    if len(payload) > MAX_REPORT_BYTES:
        raise ManifestError("ENA metadata report exceeds the 16 MiB safety limit")
    return payload


def fastq_https(location):
    """Convert a trusted ENA FTP location to HTTPS without accepting arbitrary URLs."""
    if not isinstance(location, str) or not location:
        raise ManifestError("FASTQ location is empty")
    if any(ord(c) < 33 or ord(c) == 127 for c in location) or "\\" in location:
        raise ManifestError("FASTQ location contains whitespace or control characters")
    parsed = urlsplit(location if "://" in location else "https://" + location)
    if parsed.scheme not in {"https", "ftp"} or parsed.netloc != "ftp.sra.ebi.ac.uk":
        raise ManifestError("FASTQ host must be exactly ftp.sra.ebi.ac.uk (FTP or HTTPS)")
    if parsed.query or parsed.fragment:
        raise ManifestError("FASTQ URLs must not contain queries or fragments")
    if (not re.fullmatch(r"/vol1/fastq/[A-Za-z0-9_./-]+\.f(?:ast)?q\.gz", parsed.path)
            or any(part in {"", ".", ".."} for part in parsed.path.split("/")[1:])):
        raise ManifestError("Invalid ENA FASTQ path")
    return urlunsplit(("https", "ftp.sra.ebi.ac.uk", parsed.path, "", ""))


def split_files(value, name):
    if not value:
        return []
    values = value.split(";")
    if any(not item or item != item.strip() for item in values):
        raise ManifestError(f"{name} contains empty or whitespace-padded entries")
    return values


def paper_label_mapping(record):
    """Map only concordant archive labels to this paper's explicit sample list.

    Species, tissue, individual and biological condition are deliberately absent.
    They require the paper's experimental design, not guesses from library names.
    """
    empty = {"sample_id": "", "assay": "", "mapping_evidence": "unresolved"}
    if not {record["study_accession"], record["secondary_study_accession"]} & {"PRJNA451237", "SRP142554"}:
        return empty
    name = record["library_name"]
    if name != record["experiment_alias"]:
        return empty
    if record["library_strategy"] == "RNA-Seq" and record["library_layout"] == "PAIRED" and name.startswith("mRNA_"):
        sample, assay = name[5:], "mRNA"
    elif record["library_strategy"] == "miRNA-Seq" and record["library_layout"] == "SINGLE":
        sample, assay = name, "miRNA"
        if sample == "M2-EA":
            return empty
    else:
        return empty
    if sample not in PAPER_SAMPLES:
        return empty
    return {"sample_id": sample, "assay": assay,
            "mapping_evidence": "ENA_alias_and_library_agree_with_paper_sample_label"}


def validate_record(record, expected_study):
    patterns = {
        "run_accession": r"[SED]RR\d+",
        "experiment_accession": r"[SED]RX\d+",
        "sample_accession": r"(?:SAM[END][A-Z]?\d+|[SED]RS\d+)",
        "study_accession": STUDY_PATTERN,
    }
    for field, pattern in patterns.items():
        if not re.fullmatch(pattern, record[field]):
            raise ManifestError(f"Invalid {field}: {record[field]!r}")
    if expected_study not in {record["study_accession"], record["secondary_study_accession"]}:
        raise ManifestError(f"Run {record['run_accession']} does not belong to {expected_study}")
    if record["library_layout"] not in {"SINGLE", "PAIRED"}:
        raise ManifestError("library_layout must be SINGLE or PAIRED")
    for field in ("read_count", "base_count"):
        if record[field] and not re.fullmatch(r"[0-9]+", record[field]):
            raise ManifestError(f"{field} must be a nonnegative integer or empty")
    locations = split_files(record["fastq_ftp"], "fastq_ftp")
    checksums = split_files(record["fastq_md5"], "fastq_md5")
    sizes = split_files(record["fastq_bytes"], "fastq_bytes")
    if not len(locations) == len(checksums) == len(sizes):
        raise ManifestError("fastq_ftp, fastq_md5 and fastq_bytes counts differ")
    urls = [fastq_https(location) for location in locations]
    if len(set(urls)) != len(urls):
        raise ManifestError("Duplicate FASTQ location within one run")
    if any(not re.fullmatch(r"[0-9A-Fa-f]{32}", checksum) for checksum in checksums):
        raise ManifestError("FASTQ MD5 must contain exactly 32 hexadecimal characters")
    if any(not re.fullmatch(r"[0-9]+", size) or int(size) <= 0 for size in sizes):
        raise ManifestError("FASTQ byte sizes must be positive integers")
    return {
        **record, **paper_label_mapping(record),
        "fastq_https": ";".join(urls), "fastq_file_count": str(len(urls)),
        "fastq_total_bytes": str(sum(map(int, sizes))),
    }


def parse_report(payload, study="SRP142554", retrieved_utc=None, expected_runs=None):
    source = report_url(study)
    if expected_runs is not None and expected_runs <= 0:
        raise ManifestError("expected_runs must be positive")
    if len(payload) > MAX_REPORT_BYTES:
        raise ManifestError("ENA metadata report exceeds the 16 MiB safety limit")
    try:
        raw = payload.decode("utf-8-sig")
    except UnicodeError as error:
        raise ManifestError("ENA report must be UTF-8 TSV") from error
    digest = hashlib.sha256(payload).hexdigest()
    timestamp = (datetime.now(timezone.utc).isoformat(timespec="seconds")
                 if retrieved_utc is None else retrieved_utc)
    reader = csv.DictReader(io.StringIO(raw, newline=""), delimiter="\t", strict=True)
    try:
        headers = reader.fieldnames
        if not headers or len(headers) != len(set(headers)):
            raise ManifestError("Missing or duplicate TSV column headers")
        missing = set(FIELDS) - set(headers)
        if missing:
            raise ManifestError("Missing ENA columns: " + ", ".join(sorted(missing)))
        records = []
        seen = set()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ManifestError(f"Malformed TSV row at line {reader.line_num}")
            if any(any(ord(char) < 32 or ord(char) == 127 for char in value)
                   for value in row.values()):
                raise ManifestError(f"Control character in TSV row at line {reader.line_num}")
            record = validate_record({field: row[field] for field in FIELDS}, study)
            run = record["run_accession"]
            if run in seen:
                raise ManifestError(f"Duplicate run accession: {run}")
            seen.add(run)
            records.append({**record, "source_url": source, "retrieved_utc": timestamp,
                            "source_report_sha256": digest})
    except csv.Error as error:
        raise ManifestError(f"Invalid TSV: {error}") from error
    if not records:
        raise ManifestError("ENA returned no runs; no manifest will be written")
    if expected_runs is not None and len(records) != expected_runs:
        raise ManifestError(f"Expected {expected_runs} runs, received {len(records)}")
    return sorted(records, key=lambda item: item["run_accession"])


def manifest_tsv(records):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    return stream.getvalue()


def atomic_write(output, text, overwrite=False):
    """Publish one complete TSV; failures leave any previous output untouched."""
    output = Path(output)
    if output.is_symlink():
        raise ManifestError("Refusing to replace a symlink output")
    if output.exists() and not overwrite:
        raise ManifestError(f"Output already exists: {output}; use --overwrite deliberately")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=output.parent, prefix="." + output.name + ".",
                                         suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary, output)
        else:
            # Atomic no-clobber publication, including a concurrently created output.
            os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", default="SRP142554")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input-report", type=Path, help="validate an existing ENA TSV without network access")
    parser.add_argument("--expect-runs", type=int, help="fail if the report has a different run count")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.input_report:
            with args.input_report.open("rb") as stream:
                payload = stream.read(MAX_REPORT_BYTES + 1)
        else:
            payload = fetch_report(args.study, args.timeout)
        records = parse_report(payload, args.study,
                               retrieved_utc="" if args.input_report else None,
                               expected_runs=args.expect_runs)
        atomic_write(args.output, manifest_tsv(records), args.overwrite)
    except (ManifestError, OSError, ValueError) as error:
        print(f"Metadata manifest failed: {error}", file=sys.stderr)
        return 1
    file_count = sum(int(row["fastq_file_count"]) for row in records)
    print(f"Wrote {len(records)} runs / {file_count} FASTQ references to {args.output}; no FASTQ downloaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
