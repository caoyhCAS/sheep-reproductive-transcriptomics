import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("evidence", Path(__file__).parents[1] / "scripts/filter_evidence.py")
evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evidence)


class EvidenceTests(unittest.TestCase):
    def test_ambiguous_and_policies(self):
        rows = [dict(precursor_id="a", reads="10", score="5"),
                dict(precursor_id="b", reads="11", score="4"),
                dict(precursor_id="c", reads="11", score="5"),
                dict(precursor_id="d", reads="10", score="4")]
        self.assertEqual([r["precursor_id"] for r in evidence.precursors(rows, "joint-failure")], ["a", "b", "c"])
        self.assertEqual([r["precursor_id"] for r in evidence.precursors(rows, "either-failure")], ["c"])

    def test_blast_identity_does_not_imply_coverage(self):
        rows = [dict(qseqid="q", sseqid="s", evalue="1e-10", pident="100", qcovhsp="20")]
        self.assertEqual(len(evidence.blast(rows, "identity")), 1)
        self.assertEqual(evidence.blast(rows, "coverage"), [])
        self.assertEqual(evidence.blast(rows, "both"), [])

    def test_enrichment_strict_boundary_and_correction(self):
        rows = [dict(term_id="a", adjusted_p="0.05", correction="Bonferroni"),
                dict(term_id="b", adjusted_p="0.049", correction="Bonferroni")]
        self.assertEqual([r["term_id"] for r in evidence.enrichment(rows, "Bonferroni")], ["b"])
        with self.assertRaises(ValueError):
            evidence.enrichment(rows, "BH")

    def test_invalid_values_fail(self):
        for value in ("nan", "inf", "-1", "1.01"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                evidence.enrichment([dict(term_id="a", adjusted_p=value, correction="BH")], "BH")
        with self.assertRaises(ValueError):
            evidence.precursors([dict(precursor_id="a", reads="10.5", score="5")], "either-failure")


if __name__ == "__main__":
    unittest.main()
