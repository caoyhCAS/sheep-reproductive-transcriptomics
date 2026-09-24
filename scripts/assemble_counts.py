#!/usr/bin/env python3
"""Join RSEM expected counts by ID; new reconstruction helper, not original code."""
import argparse
import csv
import math
from pathlib import Path


def load(path, columns):
    with open(path, newline="", encoding="utf-8") as handle:
        table = csv.DictReader(handle, delimiter="\t")
        if (not table.fieldnames or not set(columns) <= set(table.fieldnames)
                or len(set(table.fieldnames)) != len(table.fieldnames)):
            raise ValueError(f"{path}: invalid header; required {columns}")
        rows = list(table)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"{path}: malformed row")
    return rows


def matrix(manifest, level):
    rows = load(manifest, ["sample_id", "rsem_results"])
    if not rows:
        raise ValueError("empty manifest")
    samples, matrices, paths = [], [], set()
    for row in rows:
        sample = row["sample_id"]
        if not sample or sample != sample.strip() or sample in samples or sample == "gene_id":
            raise ValueError("invalid or duplicate sample_id")
        path = Path(row["rsem_results"])
        if not path.is_absolute():
            path = Path(manifest).resolve().parent / path
        if path.resolve() in paths:
            raise ValueError("same RSEM result used for multiple samples")
        paths.add(path.resolve())
        key = "gene_id" if level == "genes" else "transcript_id"
        values = {}
        for feature in load(path, [key, "expected_count"]):
            identifier = feature[key]
            if not identifier or identifier != identifier.strip() or identifier in values:
                raise ValueError("empty or duplicate feature ID")
            count = float(feature["expected_count"])
            if not math.isfinite(count) or count < 0:
                raise ValueError("expected counts must be finite and nonnegative")
            values[identifier] = feature["expected_count"]
        if not values:
            raise ValueError("empty RSEM result")
        if matrices and set(values) != set(matrices[0]):
            raise ValueError("feature sets differ; use one common RSEM reference")
        samples.append(sample)
        matrices.append(values)
    return samples, [[feature, *[values[feature] for values in matrices]] for feature in sorted(matrices[0])], paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="TSV sample_id,rsem_results; paths relative to manifest")
    parser.add_argument("--level", required=True, choices=["genes", "isoforms"])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        samples, rows, sources = matrix(args.manifest, args.level)
        output = Path(args.output)
        if output.exists() or output.resolve() in sources or output.resolve() == Path(args.manifest).resolve():
            raise ValueError("output must be a new file, distinct from all inputs")
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *samples])
            writer.writerows(rows)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
