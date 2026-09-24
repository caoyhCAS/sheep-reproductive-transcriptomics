#!/usr/bin/env python3
"""Method reconstruction, not original author code. Filter explicit external exports."""
import argparse
import csv
import math
from pathlib import Path


def read_table(path, required):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path}: required columns: {', '.join(required)}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("duplicate column names")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError("ragged table")
        return reader.fieldnames, rows


def number(value, label, low=None, high=None):
    result = float(value)
    if not math.isfinite(result) or (low is not None and result < low) or (high is not None and result > high):
        raise ValueError(f"invalid {label}: {value}")
    return result


def precursors(rows, policy):
    """Paper's AND wording is ambiguous: the caller must choose explicitly."""
    kept = []
    seen = set()
    for row in rows:
        identifier = row["precursor_id"]
        if not identifier or identifier in seen:
            raise ValueError("empty or duplicate precursor_id")
        seen.add(identifier)
        reads = number(row["reads"], "reads", 0)
        if not reads.is_integer():
            raise ValueError("precursor reads must be integers")
        score = number(row["score"], "score")
        low_reads, low_score = reads <= 10, score < 5
        reject = (low_reads and low_score) if policy == "joint-failure" else (low_reads or low_score)
        if not reject:
            kept.append(row)
    return kept


def blast(rows, policy):
    """BLAST export includes qcovhsp; coverage is per HSP, not union of HSPs."""
    kept = []
    for row in rows:
        if not row["qseqid"] or not row["sseqid"]:
            raise ValueError("empty BLAST identifier")
        evalue = number(row["evalue"], "evalue", 0)
        identity = number(row["pident"], "pident", 0, 100)
        coverage = number(row["qcovhsp"], "qcovhsp", 0, 100)
        matched = {"identity": identity == 100, "coverage": coverage == 100,
                   "both": identity == 100 and coverage == 100}[policy]
        if evalue <= 1e-10 and matched:
            kept.append(row)
    return kept


def enrichment(rows, correction):
    kept = []
    for row in rows:
        if not row["term_id"] or row["correction"] != correction:
            raise ValueError("term_id or exported correction does not match")
        if number(row["adjusted_p"], "adjusted_p", 0, 1) < .05:
            kept.append(row)
    return kept


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    precursor = sub.add_parser("precursors")
    precursor.add_argument("--policy", required=True, choices=["joint-failure", "either-failure"])
    blast_parser = sub.add_parser("blast")
    blast_parser.add_argument("--matching-ratio", required=True, choices=["identity", "coverage", "both"])
    enrich = sub.add_parser("enrichment")
    enrich.add_argument("--correction", required=True)
    for child in (precursor, blast_parser, enrich):
        child.add_argument("--input", required=True)
        child.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        if Path(args.input).resolve() == Path(args.output).resolve():
            raise ValueError("input and output paths must differ")
        required = {"precursors": ["precursor_id", "reads", "score"],
                    "blast": ["qseqid", "sseqid", "evalue", "pident", "qcovhsp"],
                    "enrichment": ["term_id", "adjusted_p", "correction"]}[args.command]
        header, rows = read_table(args.input, required)
        if args.command == "precursors":
            result = precursors(rows, args.policy)
        elif args.command == "blast":
            result = blast(rows, args.matching_ratio)
        else:
            result = enrichment(rows, args.correction)
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(result)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
