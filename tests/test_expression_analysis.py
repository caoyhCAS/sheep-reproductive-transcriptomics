"""Deterministic synthetic checks; standard-library unittest, also pytest compatible."""

import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "expression_analysis.py"
SPEC = importlib.util.spec_from_file_location("expression_analysis", SCRIPT)
ea = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ea
SPEC.loader.exec_module(ea)


class ExpressionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def tsv(self, name, fields, rows):
        path = self.directory / name
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(rows)
        return path

    def normalized(self, rows, *, mirna=False, name="expression.tsv"):
        return self.tsv(name, ea.NORMALIZED_FIELDS + (["class"] if mirna else []), rows)

    def rec(self, gene, fc, p=0.01, cls=None, numerator="A", denominator="B"):
        return ea.Expression(gene, fc, p, numerator, denominator, cls)

    def test_fold_change_orientation(self):
        self.assertEqual(ea.fold_change(2, 8), 2)
        self.assertEqual(ea.fold_change(8, 2), -2)

    def test_fold_change_preserves_infinity(self):
        self.assertEqual(ea.fold_change(0, 8), math.inf)
        self.assertEqual(ea.fold_change(8, 0), -math.inf)

    def test_fold_change_extreme_finite_ratio(self):
        self.assertTrue(math.isfinite(ea.fold_change(1e-300, 1e300)))

    def test_fold_change_invalid_expression(self):
        for first, second in [(0, 0), (-1, 1), (1, math.nan), (math.inf, 1)]:
            with self.subTest(pair=(first, second)), self.assertRaises(ValueError):
                ea.fold_change(first, second)

    def test_cuffdiff_selects_comparison_and_ok_status(self):
        path = self.tsv("gene_exp.diff", ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "q_value", "log2(fold_change)"],
                        [["g1", "B", "A", "OK", 1, 3.9999, 0.01, 2],
                         ["g2", "B", "A", "NOTEST", "nan", "nan", "nan", "nan"],
                         ["g1", "C", "A", "OK", 2, 16, 0.02, 3]])
        result = ea.normalize_cuffdiff(path, "A", "B")
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].log2fc, 2)
        self.assertEqual(len(ea.filter_expression(result, "mrna")), 1)
        self.assertEqual((result[0].numerator, result[0].denominator), ("A", "B"))

    def test_cuffdiff_no_automatic_reversal(self):
        path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "q_value", "log2(fold_change)"], [["g", "B", "A", "OK", 1, 8, 0.01, 3]])
        with self.assertRaisesRegex(ValueError, "no Cuffdiff rows"):
            ea.normalize_cuffdiff(path, "B", "A")

    def test_cuffdiff_duplicate_ids_rejected(self):
        path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "q_value", "log2(fold_change)"], [["g", "B", "A", "OK", 1, 8, 0.01, 3]] * 2)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            ea.normalize_cuffdiff(path, "A", "B")

    def test_mrna_and_denovo_inclusive_thresholds(self):
        records = [self.rec("yes_plus", 2), self.rec("yes_minus", -2), self.rec("low_fc", 1.99), self.rec("high_p", 3, 0.010001)]
        for mode in ("mrna", "denovo"):
            self.assertEqual([r.gene_id for r in ea.filter_expression(records, mode)], ["yes_plus", "yes_minus"])

    def test_mirna_thresholds_and_classes(self):
        records = [self.rec("a", -1, 0.05, "sheep"), self.rec("b", 1, 0.05, "conserved"), self.rec("c", 5, 0, "novel"), self.rec("d", 0.999, 0, "sheep")]
        self.assertEqual([r.gene_id for r in ea.filter_expression(records, "mirna")], ["a", "b"])

    def test_normalized_infinity_is_meaningful(self):
        path = self.normalized([["g", "-inf", 0.01, "OK", "A", "B"]])
        result = ea.filter_expression(ea.read_expression(path, "A", "B"), "mrna")
        self.assertEqual(result[0].log2fc, -math.inf)

    def test_normalized_nan_fc_rejected(self):
        path = self.normalized([["g", "nan", 0.01, "OK", "A", "B"]])
        with self.assertRaises(ValueError):
            ea.read_expression(path, "A", "B")

    def test_normalized_fc_overflow_is_not_literal_infinity(self):
        path = self.normalized([["g", "1e99999", 0.01, "OK", "A", "B"]])
        with self.assertRaisesRegex(ValueError, "overflow"):
            ea.read_expression(path, "A", "B")

    def test_cuffdiff_native_infinity_without_fpkm_columns(self):
        path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "q_value", "log2(fold_change)"], [["g", "B", "A", "OK", 0.01, "inf"]])
        self.assertEqual(ea.normalize_cuffdiff(path, "A", "B")[0].log2fc, math.inf)

    def test_cuffdiff_native_fc_required_no_fallback(self):
        path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "q_value"], [["g", "B", "A", "OK", 1, 8, 0.01]])
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            ea.normalize_cuffdiff(path, "A", "B")

    def test_cuffdiff_optional_fpkm_must_be_valid(self):
        for value in (-1, "nan", "inf"):
            path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "value_1", "q_value", "log2(fold_change)"], [["g", "B", "A", "OK", value, 0.01, 3]])
            with self.subTest(value=value), self.assertRaises(ValueError):
                ea.normalize_cuffdiff(path, "A", "B")

    def test_bad_adjusted_p_rejected(self):
        for p in ("nan", "NA", "inf", "-0.01", "1.01", "banana"):
            with self.subTest(p=p), self.assertRaises(ValueError):
                ea.adjusted_p(p)

    def test_normalized_nonok_na_skipped(self):
        path = self.normalized([["g", "NA", "NA", "NOTEST", "A", "B"]])
        self.assertEqual(ea.read_expression(path, "A", "B"), [])

    def test_unknown_status_rejected(self):
        path = self.normalized([["g", 2, 0.01, "typo", "A", "B"]])
        with self.assertRaisesRegex(ValueError, "unknown status"):
            ea.read_expression(path, "A", "B")

    def test_normalized_duplicate_id_rejected(self):
        path = self.normalized([["g", 2, 0.01, "OK", "A", "B"]] * 2)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            ea.read_expression(path, "A", "B")

    def test_normalized_orientation_mismatch_rejected(self):
        path = self.normalized([["g", 2, 0.01, "OK", "B", "A"]])
        with self.assertRaisesRegex(ValueError, "orientation"):
            ea.read_expression(path, "A", "B")

    def test_nonok_orientation_mismatch_still_rejected(self):
        path = self.normalized([["g", "NA", "NA", "NOTEST", "B", "A"]])
        with self.assertRaisesRegex(ValueError, "orientation"):
            ea.read_expression(path, "A", "B")

    def test_mirna_class_required_and_checked(self):
        path = self.normalized([["g", 2, 0.01, "OK", "A", "B"]])
        with self.assertRaisesRegex(ValueError, "class"):
            ea.read_expression(path, "A", "B", mirna=True)
        path = self.normalized([["g", 2, 0.01, "OK", "A", "B", "unknown"]], mirna=True)
        with self.assertRaisesRegex(ValueError, "class"):
            ea.read_expression(path, "A", "B", mirna=True)

    def test_strict_headers_and_row_width(self):
        for text in ("gene_id\tgene_id\ng\tg\n", "gene_id\tpadj\ng\n", "gene_id\tpadj\ng\t0.1\textra\n", " gene_id\np\n"):
            path = self.directory / "bad.tsv"
            path.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                ea.read_tsv(path, ["gene_id"])

    def test_target_dedup_uses_maximum_score(self):
        path = self.tsv("targets.tsv", ["mirna", "gene_id", "context_percentile"], [["mi", "g", 40], ["mi", "g", 60], ["mi", "g", 50]])
        self.assertEqual(ea.read_targets(path), {("mi", "g"): 60})

    def test_target_percentile_invalid(self):
        for value in (-1, 101, "nan", "inf", "invalid"):
            path = self.tsv("targets.tsv", ["mirna", "gene_id", "context_percentile"], [["mi", "g", value]])
            with self.subTest(value=value), self.assertRaises(ValueError):
                ea.read_targets(path)

    def test_network_requires_opposite_signs_and_threshold(self):
        genes = [self.rec("up", 2), self.rec("down", -2), self.rec("zero", 0)]
        mirnas = [self.rec("mi_down", -1, 0.05, "sheep"), self.rec("mi_up", 1, 0.05, "conserved")]
        targets = {("mi_down", "up"): 50, ("mi_down", "down"): 80, ("mi_up", "down"): 80,
                   ("mi_up", "up"): 90, ("mi_down", "zero"): 90}
        network, core = ea.integrate(genes, mirnas, targets, "A", "B")
        self.assertEqual([(r["mirna"], r["gene_id"]) for r in network], [("mi_down", "up"), ("mi_up", "down")])
        self.assertEqual(len(core), 1)
        self.assertEqual(core[0]["core_mirna_rank"], 1)

    def test_network_reapplies_de_and_class_filters(self):
        genes = [self.rec("g", 3, 0.01), self.rec("nonde", 1, 0)]
        mirnas = [self.rec("novel", -3, 0, "novel"), self.rec("nonde", -0.5, 0, "sheep"), self.rec("ok", -3, 0, "conserved")]
        network, _ = ea.integrate(genes, mirnas, {(mi.gene_id, gene.gene_id): 80 for mi in mirnas for gene in genes}, "A", "B")
        self.assertEqual([(r["mirna"], r["gene_id"]) for r in network], [("ok", "g")])

    def test_network_rejects_mismatched_orientation(self):
        with self.assertRaisesRegex(ValueError, "orientation"):
            ea.integrate([self.rec("g", 2)], [self.rec("mi", -2, cls="sheep", numerator="B", denominator="A")], {("mi", "g"): 90}, "A", "B")

    def test_core_top10_rank_and_stable_ties(self):
        genes = [self.rec("up", 2)]
        mirnas = [self.rec(f"mi{i:02}", -(i + 1), cls="sheep") for i in range(12)]
        mirnas += [self.rec("tieA", -20, cls="conserved"), self.rec("tieB", -20, cls="conserved"), self.rec("no_target", -100, cls="sheep")]
        targets = {(mi.gene_id, "up"): 90 for mi in mirnas if mi.gene_id != "no_target"}
        network, core = ea.integrate(genes, mirnas, targets, "A", "B")
        self.assertEqual(len(network), 14)
        self.assertEqual(len(core), 10)
        self.assertEqual([r["mirna"] for r in core[:3]], ["tieA", "tieB", "mi11"])
        self.assertEqual([r["core_mirna_rank"] for r in core], list(range(1, 11)))

    def test_infinite_opposite_fold_changes_network(self):
        network, core = ea.integrate([self.rec("g", math.inf)], [self.rec("mi", -math.inf, cls="sheep")], {("mi", "g"): 50}, "A", "B")
        self.assertEqual(len(network), 1)
        self.assertEqual(len(core), 1)

    def test_empty_network_is_valid(self):
        self.assertEqual(ea.integrate([], [], {}, "A", "B"), ([], []))

    def test_empty_and_equal_orientation_labels_rejected(self):
        for pair in [("A", "A"), ("", "B"), ("A ", "B")]:
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                ea.validate_orientation(*pair)

    def test_input_output_and_duplicate_output_paths_rejected(self):
        with self.assertRaises(ValueError):
            ea.validate_paths([str(self.directory / "input")], [str(self.directory / "input")])
        with self.assertRaises(ValueError):
            ea.validate_paths([], [str(self.directory / "out"), str(self.directory / "out")])

    def test_cli_end_to_end_synthetic_network(self):
        cuffdiff = self.tsv("gene_exp.diff", ["gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "q_value", "log2(fold_change)"], [["g1", "B", "A", "OK", 1, 8, 0.005, 3]])
        mrna = self.directory / "normalized" / "mrna.tsv"
        self.assertEqual(ea.main(["normalize-cuffdiff", "--input", str(cuffdiff), "--output", str(mrna), "--numerator", "A", "--denominator", "B"]), 0)
        mirna = self.normalized([["mi1", -1, 0.05, "OK", "A", "B", "sheep"]], mirna=True, name="mirna.tsv")
        targets = self.tsv("targets.tsv", ["mirna", "gene_id", "context_percentile"], [["mi1", "g1", 50]])
        output, core = self.directory / "network.tsv", self.directory / "core.tsv"
        self.assertEqual(ea.main(["integrate", "--mrna", str(mrna), "--mirna", str(mirna), "--targets", str(targets), "--output", str(output), "--core-output", str(core), "--numerator", "A", "--denominator", "B"]), 0)
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["mrna_log2fc"], "3.0")
        self.assertEqual(rows[0]["core_mirna_rank"], "1")

    def test_cli_requires_orientation_flags(self):
        completed = subprocess.run([sys.executable, str(SCRIPT), "filter", "--input", "missing", "--output", "unused", "--mode", "mrna"], capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--numerator", completed.stderr)

    def test_filter_preserves_nonindependent_library_scope(self):
        scope = "descriptive_nonindependent_libraries"
        path = self.tsv("scoped.tsv", ea.NORMALIZED_FIELDS + ["analysis_scope"],
                        [["g", 3, 0.005, "OK", "A", "B", scope]])
        output = self.directory / "filtered.tsv"
        ea.main(["filter", "--input", str(path), "--output", str(output), "--mode", "mrna",
                 "--numerator", "A", "--denominator", "B"])
        rows = ea.read_tsv(output, ["analysis_scope"])
        self.assertEqual(rows[0]["analysis_scope"], scope)
        self.assertEqual(ea.read_expression(output, "A", "B")[0].analysis_scope, scope)

    def test_missing_analysis_scope_explicitly_unspecified(self):
        path = self.normalized([["g", 3, 0.005, "OK", "A", "B"]])
        records = ea.read_expression(path, "A", "B")
        self.assertEqual(records[0].analysis_scope, "unspecified")
        self.assertEqual(records[0].as_row()["analysis_scope"], "unspecified")

    def test_network_and_core_keep_separate_analysis_scopes(self):
        mrna = [ea.Expression("g", 3, 0.005, "A", "B", analysis_scope="descriptive_nonindependent_libraries")]
        mirna = [ea.Expression("mi", -2, 0.01, "A", "B", "sheep", "unspecified")]
        network, core = ea.integrate(mrna, mirna, {("mi", "g"): 90}, "A", "B")
        for row in (network[0], core[0]):
            self.assertEqual(row["mrna_analysis_scope"], "descriptive_nonindependent_libraries")
            self.assertEqual(row["mirna_analysis_scope"], "unspecified")

    def test_cuffdiff_scope_is_opt_in_and_default_unspecified(self):
        path = self.tsv("diff.tsv", ["gene_id", "sample_1", "sample_2", "status", "q_value", "log2(fold_change)"],
                        [["g", "B", "A", "OK", 0.005, 3]])
        self.assertEqual(ea.normalize_cuffdiff(path, "A", "B")[0].analysis_scope, "unspecified")
        output = self.directory / "normalized.tsv"
        ea.main(["normalize-cuffdiff", "--input", str(path), "--output", str(output),
                 "--numerator", "A", "--denominator", "B",
                 "--analysis-scope", "descriptive_nonindependent_libraries"])
        self.assertEqual(ea.read_tsv(output, ["analysis_scope"])[0]["analysis_scope"], "descriptive_nonindependent_libraries")


if __name__ == "__main__":
    unittest.main()
