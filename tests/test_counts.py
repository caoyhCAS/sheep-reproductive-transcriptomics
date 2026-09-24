import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("counts", Path(__file__).parents[1] / "scripts/assemble_counts.py")
counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(counts)


class CountTests(unittest.TestCase):
    def fixture(self, directory, second):
        root = Path(directory)
        (root / "a.tsv").write_text("gene_id\texpected_count\ng2\t0\ng1\t1.5\n")
        (root / "b.tsv").write_text(second)
        (root / "manifest.tsv").write_text("sample_id\trsem_results\nS1\ta.tsv\nS2\tb.tsv\n")
        return root / "manifest.tsv"

    def test_ids_aligned_and_fractional_counts_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = self.fixture(directory, "gene_id\texpected_count\ng1\t2.25\ng2\t3\n")
            samples, rows, _ = counts.matrix(manifest, "genes")
            self.assertEqual(samples, ["S1", "S2"])
            self.assertEqual(rows, [["g1", "1.5", "2.25"], ["g2", "0", "3"]])

    def test_mismatched_or_duplicate_features_rejected(self):
        for data in ["gene_id\texpected_count\ng1\t2\n", "gene_id\texpected_count\ng1\t2\ng1\t3\n",
                     "gene_id\texpected_count\ng1\tnan\ng2\t3\n"]:
            with self.subTest(data=data), tempfile.TemporaryDirectory() as directory:
                manifest = self.fixture(directory, data)
                with self.assertRaises(ValueError):
                    counts.matrix(manifest, "genes")
