"""Base-R input contracts only; never install or mock statistical packages.

Tests skip explicitly if Rscript is absent. They do not test edgeR/DESeq numerics.
"""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript unavailable: R input contracts were not executed")
class RInputContracts(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="r contracts ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.counts = self.directory / "counts.tsv"
        self.samples = self.directory / "samples.tsv"
        self.features = self.directory / "features.tsv"
        self.output = self.directory / "new results" / "result.tsv"
        self.counts.write_text("gene_id\ts1\ts2\ts3\ts4\ng1\t10\t20\t30\t40\ng2\t25\t20\t8\t6\n")
        self.samples.write_text("sample_id\tcondition\tindividual_id\ns1\tden\ti1\ns2\tden\ti2\ns3\tnum\ti3\ns4\tnum\ti4\n")
        self.features.write_text("gene_id\tclass\ng1\tsheep\ng2\tconserved\n")

    def command(self, script, additional=(), validate=True):
        args = [
            RSCRIPT, str(ROOT / "scripts" / script),
            "--counts", str(self.counts), "--samples", str(self.samples),
            "--numerator", "num", "--denominator", "den", "--output", str(self.output),
        ]
        if script == "deseq_mirna.R":
            args += ["--features", str(self.features), "--dispersion-method", "pooled", "--sharing-mode", "maximum", "--fit-type", "parametric"]
        if validate:
            args.append("--validate-only")
        return args + list(additional)

    def run_script(self, script, additional=(), validate=True):
        return subprocess.run(self.command(script, additional, validate), text=True, capture_output=True, check=False)

    def assert_rejected(self, result, message):
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(message, result.stderr)
        self.assertFalse(self.output.exists())

    def test_both_scripts_parse_and_validate_without_packages_or_outputs(self):
        for script in ("edger_denovo.R", "deseq_mirna.R"):
            with self.subTest(script=script):
                result = self.run_script(script)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("VALIDATION ONLY", result.stdout)
                self.assertFalse(self.output.parent.exists())

    def test_fractional_rsem_counts_accepted_only_by_edger(self):
        self.counts.write_text(self.counts.read_text().replace("\t10\t", "\t10.25\t"))
        result = self.run_script("edger_denovo.R")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_rejected(self.run_script("deseq_mirna.R"), "integer read counts")

    def test_m2_repeated_material_requires_explicit_descriptive_override(self):
        self.samples.write_text("sample_id\tcondition\tindividual_id\ns1\tden\tM2\ns2\tden\tM2\ns3\tnum\tM2\ns4\tnum\tM2\n")
        for script in ("edger_denovo.R", "deseq_mirna.R"):
            with self.subTest(script=script):
                self.assert_rejected(self.run_script(script), "--allow-pseudoreplication")
                result = self.run_script(script, ["--allow-pseudoreplication"])
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("DESCRIPTIVE ONLY", result.stderr)

    def test_padding_cannot_turn_one_animal_into_two(self):
        self.samples.write_text("sample_id\tcondition\tindividual_id\ns1\tden\tM2\ns2\tden\t M2\ns3\tnum\tM2\ns4\tnum\tM2 \n")
        self.assert_rejected(self.run_script("edger_denovo.R"), "individual_id must not contain padding")

    def test_nonfinite_and_negative_counts_rejected(self):
        original = self.counts.read_text()
        for bad in ("-1", "NaN", "Inf", "NA"):
            with self.subTest(value=bad):
                self.counts.write_text(original.replace("\t10\t", "\t" + bad + "\t"))
                self.assert_rejected(self.run_script("edger_denovo.R"), "finite, nonnegative")

    def test_duplicate_features_rejected(self):
        self.counts.write_text(self.counts.read_text().replace("g2\t", "g1\t"))
        self.assert_rejected(self.run_script("edger_denovo.R"), "gene_id values must be nonempty and unique")

    def test_sample_mismatch_is_not_silently_ignored(self):
        self.samples.write_text(self.samples.read_text().replace("s4\t", "absent\t"))
        self.assert_rejected(self.run_script("deseq_mirna.R"), "must match exactly")

    def test_metadata_can_be_reordered(self):
        lines = self.samples.read_text().splitlines()
        self.samples.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")
        result = self.run_script("edger_denovo.R")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_numerator_condition_rejected(self):
        self.samples.write_text(self.samples.read_text().replace("\tnum\t", "\tother\t"))
        self.assert_rejected(self.run_script("edger_denovo.R"), "Both numerator and denominator")

    def test_missing_individual_column_rejected(self):
        self.samples.write_text(self.samples.read_text().replace("individual_id", "technical_library"))
        self.assert_rejected(self.run_script("edger_denovo.R"), "sample_id, condition, individual_id")

    def test_no_replication_requires_explicit_dispersion_choice(self):
        self.counts.write_text("gene_id\ts1\ts3\ng1\t10\t30\ng2\t25\t8\n")
        self.samples.write_text("sample_id\tcondition\tindividual_id\ns1\tden\ti1\ns3\tnum\ti3\n")
        self.assert_rejected(self.run_script("edger_denovo.R", ["--allow-pseudoreplication"]), "--dispersion")
        result = self.run_script("edger_denovo.R", ["--allow-pseudoreplication", "--dispersion", "0.2"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_rejected(self.run_script("deseq_mirna.R", ["--allow-pseudoreplication"]), "blind dispersion")
        args = self.command("deseq_mirna.R", ["--allow-pseudoreplication"])
        args[args.index("pooled")] = "blind"
        result = subprocess.run(args, text=True, capture_output=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_mirna_feature_classes_and_ids_validated(self):
        self.features.write_text("gene_id\tclass\ng1\tunknown\ng2\tconserved\n")
        self.assert_rejected(self.run_script("deseq_mirna.R"), "class must be exactly")
        self.features.write_text("gene_id\tclass\ng1\tsheep\nother\tconserved\n")
        self.assert_rejected(self.run_script("deseq_mirna.R"), "Feature IDs must exactly match")

    def test_novel_only_features_are_not_tested(self):
        self.features.write_text("gene_id\tclass\ng1\tnovel\ng2\tnovel\n")
        self.assert_rejected(self.run_script("deseq_mirna.R"), "No nonzero sheep/conserved")

    def test_no_implicit_pseudocount_for_size_factors(self):
        self.counts.write_text("gene_id\ts1\ts2\ts3\ts4\ng1\t10\t20\t0\t0\ng2\t0\t0\t8\t6\n")
        self.assert_rejected(self.run_script("deseq_mirna.R"), "no pseudocount is added")

    def test_existing_output_is_not_modified(self):
        self.output.parent.mkdir()
        self.output.write_text("existing user result\n")
        result = self.run_script("edger_denovo.R")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Refusing to overwrite", result.stderr)
        self.assertEqual(self.output.read_text(), "existing user result\n")

    def test_unknown_and_duplicate_flags_rejected(self):
        self.assert_rejected(self.run_script("edger_denovo.R", ["--unreported-setting", "1"]), "Unknown option")
        self.assert_rejected(self.run_script("edger_denovo.R", ["--numerator", "reverse"]), "Duplicate option")


if __name__ == "__main__":
    unittest.main()
