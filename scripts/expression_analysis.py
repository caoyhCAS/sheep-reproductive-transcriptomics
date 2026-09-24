#!/usr/bin/env python3
"""Method reconstruction for DOI 10.1038/s41437-018-0090-1; not original author code.

Normalize Cuffdiff, filter exported differential-expression results, and integrate
exported TargetScan 7 scores. This script does not run target prediction or fit a
differential-expression model. All log2 fold changes use numerator/denominator.

Reconstruction decisions: Cuffdiff's native log2(fold_change) is preserved with
sample_2/sample_1 direction; non-OK rows are excluded; duplicate IDs are rejected;
the core network uses
the ten most downregulated *network-participating* eligible miRNAs, sorted by
increasing log2FC then ID. The paper does not specify this ranking/tie convention.
Infinite native log2FC is retained; NaN and invalid adjusted P values are rejected
for OK rows. Native fold changes are never recomputed from rounded FPKM values.
Analysis-scope labels are carried through filtering and separately for each RNA
type in network outputs. Missing scope is explicitly "unspecified"; passing a
statistical threshold never establishes validity of biological replication.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys


NORMALIZED_FIELDS = ["gene_id", "log2fc", "padj", "status", "numerator", "denominator"]
NORMALIZED_OUTPUT_FIELDS = NORMALIZED_FIELDS + ["analysis_scope"]
NETWORK_FIELDS = ["mirna", "gene_id", "context_percentile", "mirna_log2fc", "mirna_padj",
                  "mrna_log2fc", "mrna_padj", "mirna_class", "numerator", "denominator",
                  "core_mirna_rank", "mrna_analysis_scope", "mirna_analysis_scope"]
VALID_STATUSES = {"OK", "NOTEST", "LOWDATA", "HIDATA", "FAIL"}
VALID_CLASSES = {"sheep", "conserved", "novel"}
THRESHOLDS = {"mrna": (2.0, 0.01), "denovo": (2.0, 0.01), "mirna": (1.0, 0.05)}


@dataclass(frozen=True)
class Expression:
    gene_id: str
    log2fc: float
    padj: float
    numerator: str
    denominator: str
    mirna_class: str | None = None
    analysis_scope: str = "unspecified"

    def as_row(self, *, with_class: bool = False) -> dict:
        row = {"gene_id": self.gene_id, "log2fc": self.log2fc, "padj": self.padj,
               "status": "OK", "numerator": self.numerator, "denominator": self.denominator,
               "analysis_scope": checked_scope(self.analysis_scope)}
        if with_class:
            row["class"] = self.mirna_class
        return row


def identifier(value: str, field: str) -> str:
    if not value or value != value.strip():
        raise ValueError(f"{field} must be nonempty and have no surrounding whitespace")
    return value


def validate_orientation(numerator: str, denominator: str) -> None:
    identifier(numerator, "numerator")
    identifier(denominator, "denominator")
    if numerator == denominator:
        raise ValueError("numerator and denominator must differ")


def numeric(value: str, field: str, *, allow_infinity: bool = False) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{field}: invalid numeric value {value!r}") from error
    if math.isnan(number) or (not allow_infinity and not math.isfinite(number)):
        raise ValueError(f"{field}: NaN/nonfinite value is not permitted")
    if math.isinf(number) and value.strip().lower() not in {"inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}:
        raise ValueError(f"{field}: numeric overflow is not an explicit infinite statistic")
    return number


def adjusted_p(value: str) -> float:
    result = numeric(value, "padj")
    if not 0 <= result <= 1:
        raise ValueError("padj must lie in [0,1]")
    return result


def read_tsv(path: str | Path, required: list[str]) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)) or any(not key or key != key.strip() for key in header):
            raise ValueError(f"{path}: missing, duplicate, empty, or whitespace-padded header")
        absent = set(required).difference(header)
        if absent:
            raise ValueError(f"{path}: missing required columns: {', '.join(sorted(absent))}")
        rows = []
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}, line {line}: wrong number of TSV fields")
            rows.append(row)
    return rows


def write_tsv(path: str | Path, fields: list[str], rows: list[dict]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fold_change(value_1: float, value_2: float) -> float:
    """log2(value_2/value_1), preserving +/-infinity but rejecting undefined 0/0."""
    if any(not math.isfinite(value) or value < 0 for value in (value_1, value_2)):
        raise ValueError("Cuffdiff expression values must be finite and nonnegative")
    if value_1 == value_2 == 0:
        raise ValueError("undefined Cuffdiff fold change: value_1 and value_2 are both zero")
    if value_1 == 0:
        return math.inf
    if value_2 == 0:
        return -math.inf
    # Subtract logs to avoid overflow/underflow in the expression ratio itself.
    return math.log2(value_2) - math.log2(value_1)


def checked_status(value: str) -> str:
    if value not in VALID_STATUSES:
        raise ValueError(f"unknown status {value!r}; expected one of {', '.join(sorted(VALID_STATUSES))}")
    return value


def checked_scope(value: str | None) -> str:
    """Unknown scope remains unknown; this is provenance, not a validity verdict."""
    if value is None or not value.strip():
        return "unspecified"
    return identifier(value, "analysis_scope")


def normalize_cuffdiff(path: str | Path, numerator: str, denominator: str,
                      analysis_scope: str = "unspecified") -> list[Expression]:
    validate_orientation(numerator, denominator)
    scope = checked_scope(analysis_scope)
    rows = read_tsv(path, ["gene_id", "sample_1", "sample_2", "status", "log2(fold_change)", "q_value"])
    output = []
    seen = set()
    selected = 0
    for row in rows:
        # A Cuffdiff file may hold many comparisons; select the exact requested one.
        if row["sample_1"] != denominator or row["sample_2"] != numerator:
            continue
        selected += 1
        gene = identifier(row["gene_id"], "gene_id")
        if gene in seen:
            raise ValueError(f"duplicate gene_id in requested Cuffdiff comparison: {gene}")
        seen.add(gene)
        if checked_status(row["status"]) != "OK":
            continue
        for expression_field in ("value_1", "value_2"):
            if expression_field in row and numeric(row[expression_field], expression_field) < 0:
                raise ValueError("Cuffdiff expression values must be nonnegative")
        output.append(Expression(
            gene,
            numeric(row["log2(fold_change)"], "log2(fold_change)", allow_infinity=True),
            adjusted_p(row["q_value"]), numerator, denominator, analysis_scope=scope,
        ))
    if not selected:
        raise ValueError("no Cuffdiff rows match sample_2=numerator and sample_1=denominator; orientation is not reversed automatically")
    return output


def read_expression(path: str | Path, numerator: str, denominator: str, *, mirna: bool = False) -> list[Expression]:
    validate_orientation(numerator, denominator)
    rows = read_tsv(path, NORMALIZED_FIELDS + (["class"] if mirna else []))
    result = []
    seen = set()
    for row in rows:
        gene = identifier(row["gene_id"], "gene_id")
        if gene in seen:
            raise ValueError(f"duplicate gene_id in normalized table: {gene}")
        seen.add(gene)
        if (row["numerator"], row["denominator"]) != (numerator, denominator):
            raise ValueError(f"{gene}: comparison orientation does not match numerator/denominator arguments")
        mirna_class = row.get("class") if mirna else None
        if mirna and mirna_class not in VALID_CLASSES:
            raise ValueError(f"{gene}: miRNA class must be sheep, conserved, or novel")
        if checked_status(row["status"]) != "OK":
            continue
        result.append(Expression(gene, numeric(row["log2fc"], "log2fc", allow_infinity=True),
                                 adjusted_p(row["padj"]), numerator, denominator, mirna_class,
                                 checked_scope(row.get("analysis_scope"))))
    return result


def filter_expression(records: list[Expression], mode: str) -> list[Expression]:
    if mode not in THRESHOLDS:
        raise ValueError("mode must be mrna, denovo, or mirna")
    minimum_fc, maximum_p = THRESHOLDS[mode]
    result = []
    for record in records:
        # Validate public library inputs as strictly as parsed input.
        if math.isnan(record.log2fc) or not math.isfinite(record.padj) or not 0 <= record.padj <= 1:
            raise ValueError("expression record has NaN fold change or invalid adjusted P")
        if mode == "mirna":
            if record.mirna_class not in VALID_CLASSES:
                raise ValueError("miRNA class must be sheep, conserved, or novel")
            if record.mirna_class == "novel":
                continue
        if abs(record.log2fc) >= minimum_fc and record.padj <= maximum_p:
            result.append(record)
    return result


def read_targets(path: str | Path) -> dict[tuple[str, str], float]:
    """Consume normalized, already exported TargetScan scores; never predict targets."""
    rows = read_tsv(path, ["mirna", "gene_id", "context_percentile"])
    pairs = {}
    for row in rows:
        mirna = identifier(row["mirna"], "mirna")
        gene = identifier(row["gene_id"], "gene_id")
        percentile = numeric(row["context_percentile"], "context_percentile")
        if not 0 <= percentile <= 100:
            raise ValueError("context_percentile must lie in [0,100], not a 0-to-1 fractional score")
        key = (mirna, gene)
        pairs[key] = max(pairs.get(key, -math.inf), percentile)
    return pairs


def integrate(mrna: list[Expression], mirna: list[Expression], targets: dict[tuple[str, str], float],
              numerator: str, denominator: str) -> tuple[list[dict], list[dict]]:
    validate_orientation(numerator, denominator)
    # Recheck direction and duplicates for callers using the library directly.
    for records in (mrna, mirna):
        seen = set()
        for record in records:
            if record.gene_id in seen:
                raise ValueError(f"duplicate expression identifier: {record.gene_id}")
            seen.add(record.gene_id)
            if (record.numerator, record.denominator) != (numerator, denominator):
                raise ValueError("mRNA and miRNA tables must use the exact same explicit comparison orientation")
    mrna_de = {record.gene_id: record for record in filter_expression(mrna, "mrna")}
    mirna_de = {record.gene_id: record for record in filter_expression(mirna, "mirna")}
    network = []
    for (mirna_id, gene_id), percentile in sorted(targets.items()):
        if not math.isfinite(percentile) or not 0 <= percentile <= 100:
            raise ValueError("context_percentile must be finite and lie in [0,100]")
        if percentile < 50 or mirna_id not in mirna_de or gene_id not in mrna_de:
            continue
        mi, gene = mirna_de[mirna_id], mrna_de[gene_id]
        # Explicit signs avoid inf*0 -> NaN and make zero exclusion clear.
        inverse = (mi.log2fc < 0 < gene.log2fc) or (gene.log2fc < 0 < mi.log2fc)
        if not inverse:
            continue
        network.append({
            "mirna": mirna_id, "gene_id": gene_id, "context_percentile": percentile,
            "mirna_log2fc": mi.log2fc, "mirna_padj": mi.padj,
            "mrna_log2fc": gene.log2fc, "mrna_padj": gene.padj,
            "mirna_class": mi.mirna_class, "numerator": numerator, "denominator": denominator,
            "core_mirna_rank": "NA",
            "mrna_analysis_scope": checked_scope(gene.analysis_scope),
            "mirna_analysis_scope": checked_scope(mi.analysis_scope),
        })
    participating_down = {row["mirna"] for row in network if row["mirna_log2fc"] < 0}
    ranked = sorted(participating_down, key=lambda key: (mirna_de[key].log2fc, key))[:10]
    ranks = {key: rank for rank, key in enumerate(ranked, start=1)}
    for row in network:
        row["core_mirna_rank"] = ranks.get(row["mirna"], "NA")
    core = [row.copy() for row in network if row["mirna"] in ranks]
    core.sort(key=lambda row: (row["core_mirna_rank"], row["gene_id"]))
    return network, core


def validate_paths(inputs: list[str], outputs: list[str]) -> None:
    resolved_inputs = {Path(path).resolve() for path in inputs}
    resolved_outputs = [Path(path).resolve() for path in outputs]
    if len(resolved_outputs) != len(set(resolved_outputs)):
        raise ValueError("output paths must differ from each other")
    if resolved_inputs.intersection(resolved_outputs):
        raise ValueError("input and output paths must differ")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    normalize = commands.add_parser("normalize-cuffdiff", help="Select exact Cuffdiff comparison and normalize OK rows")
    normalize.add_argument("--input", required=True)
    normalize.add_argument("--output", required=True)
    normalize.add_argument("--analysis-scope", default="unspecified",
                           help="Carry supplied design/provenance label into outputs; default unspecified")
    filtering = commands.add_parser("filter", help="Apply paper's inclusive DE thresholds and eligible miRNA classes")
    filtering.add_argument("--input", required=True)
    filtering.add_argument("--output", required=True)
    filtering.add_argument("--mode", required=True, choices=sorted(THRESHOLDS))
    integration = commands.add_parser("integrate", help="Integrate exported TargetScan scores with oppositely regulated DE pairs")
    integration.add_argument("--mrna", required=True)
    integration.add_argument("--mirna", required=True)
    integration.add_argument("--targets", required=True)
    integration.add_argument("--output", required=True)
    integration.add_argument("--core-output", required=True)
    for command in (normalize, filtering, integration):
        command.add_argument("--numerator", required=True, help="Exact numerator label of every log2FC")
        command.add_argument("--denominator", required=True, help="Exact denominator label of every log2FC")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_orientation(args.numerator, args.denominator)
        if args.command == "integrate":
            validate_paths([args.mrna, args.mirna, args.targets], [args.output, args.core_output])
            mrna = read_expression(args.mrna, args.numerator, args.denominator)
            mirna = read_expression(args.mirna, args.numerator, args.denominator, mirna=True)
            network, core = integrate(mrna, mirna, read_targets(args.targets), args.numerator, args.denominator)
            write_tsv(args.output, NETWORK_FIELDS, network)
            write_tsv(args.core_output, NETWORK_FIELDS, core)
        else:
            validate_paths([args.input], [args.output])
            if args.command == "normalize-cuffdiff":
                records = normalize_cuffdiff(args.input, args.numerator, args.denominator, args.analysis_scope)
                with_class = False
            else:
                with_class = args.mode == "mirna"
                records = filter_expression(read_expression(args.input, args.numerator, args.denominator,
                                                            mirna=with_class), args.mode)
            fields = NORMALIZED_OUTPUT_FIELDS + (["class"] if with_class else [])
            write_tsv(args.output, fields, [record.as_row(with_class=with_class) for record in records])
    except (ValueError, OSError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
