#!/usr/bin/env python3
"""Offline synthetic contract demonstration. Contains no published study results."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


def write(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/demo")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("expression_analysis.py")
    direction = ["--numerator", "Finnsheep", "--denominator", "Mouflon"]
    write(output / "cuffdiff.tsv",
          ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "log2(fold_change)", "q_value"],
          [["synthetic_gene_up", "Mouflon", "Finnsheep", "OK", 7, 28, 2, .01],
           ["synthetic_gene_down", "Mouflon", "Finnsheep", "OK", 32, 4, -3, .001],
           ["synthetic_gene_ns", "Mouflon", "Finnsheep", "OK", 1, 2, 1, .2]])
    write(output / "mirna.tsv", ["gene_id", "log2fc", "padj", "status", "numerator", "denominator", "class"],
          [["synthetic_mi_down", -1, .05, "OK", "Finnsheep", "Mouflon", "sheep"],
           ["synthetic_mi_up", 2, .001, "OK", "Finnsheep", "Mouflon", "conserved"],
           ["synthetic_mi_novel", -3, .001, "OK", "Finnsheep", "Mouflon", "novel"]])
    write(output / "targets.tsv", ["mirna", "gene_id", "context_percentile"],
          [["synthetic_mi_down", "synthetic_gene_up", 50],
           ["synthetic_mi_up", "synthetic_gene_down", 90],
           ["synthetic_mi_novel", "synthetic_gene_up", 99],
           ["synthetic_mi_up", "synthetic_gene_up", 99]])
    commands = [
        ["normalize-cuffdiff", "--input", str(output / "cuffdiff.tsv"), "--output", str(output / "mrna.tsv")],
        ["filter", "--input", str(output / "mrna.tsv"), "--output", str(output / "mrna_de.tsv"), "--mode", "mrna"],
        ["filter", "--input", str(output / "mirna.tsv"), "--output", str(output / "mirna_de.tsv"), "--mode", "mirna"],
        ["integrate", "--mrna", str(output / "mrna_de.tsv"), "--mirna", str(output / "mirna_de.tsv"),
         "--targets", str(output / "targets.tsv"), "--output", str(output / "network.tsv"), "--core-output", str(output / "core.tsv")]]
    for command in commands:
        subprocess.run([sys.executable, str(script), *command, *direction], check=True)
    with (output / "network.tsv").open() as handle:
        network = list(csv.DictReader(handle, delimiter="\t"))
    with (output / "core.tsv").open() as handle:
        core = list(csv.DictReader(handle, delimiter="\t"))
    pairs = {(row["mirna"], row["gene_id"]) for row in network}
    expected = {("synthetic_mi_down", "synthetic_gene_up"), ("synthetic_mi_up", "synthetic_gene_down")}
    if pairs != expected or len(core) != 1 or core[0]["mirna"] != "synthetic_mi_down":
        raise RuntimeError("synthetic integration result differs from the expected two inverse edges")
    summary = {"data": "synthetic, not study results", "network_edges": len(network),
               "core_edges": len(core), "validation": "passed", "numerator": "Finnsheep", "denominator": "Mouflon"}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
