# Validation record

Local verification date: 2026-09-24. Version: 0.1.0.

- Python test suite: **129 passed**, **16 skipped**. The skipped tests require
  `Rscript`, which is absent from this runtime; they are not counted as passes.
- Synthetic command-line demo: passed, yielding the expected two inverse
  miRNA–mRNA edges and one core-network edge.
- Default all-stage plan: generated 126 command steps and 26 explicit
  unresolved-configuration/acknowledgment blockers; no external biology tool
  was invoked. These blockers are intended configuration requirements.
- Source Python files compiled successfully.
- Independent review checked sample/individual identities, native Cuffdiff
  fold-change preservation, inclusive DE thresholds, TargetScan percentile
  boundary, contrast directions, input alias checks, overwrite protection,
  fractional RSEM counts and analysis-scope propagation.
- Verified ENA metadata snapshot: 31 runs and 47 FASTQ references; no raw
  FASTQ file was downloaded or checksum-verified locally.

The local test environment used Python 3.12.14 and pytest 8.3.5. GitHub CI is
configured for Python 3.11 and base R 4.4.3 so the R input-contract tests can run
there. This base R installation is for validation only; it is not the paper's
historical statistical environment.

The source is published at
[caoyhCAS/sheep-reproductive-transcriptomics](https://github.com/caoyhCAS/sheep-reproductive-transcriptomics).
The local counts above record pre-publication validation. The hosted workflow
also executes the 16 base-R input-contract tests. Consult the relevant commit's
[GitHub Actions result](https://github.com/caoyhCAS/sheep-reproductive-transcriptomics/actions/workflows/ci.yml)
for its actual pass/fail counts; a configured test is not itself a passing test.
Deployment review corrected the R output-path guard to preserve `NA` from
`Sys.readlink` for nonexistent files, preventing fresh output paths from being
incorrectly rejected while keeping existing-file/symlink protection.

No TopHat/Cufflinks/Trinity/RSEM/miRDeep2 biological computation, edgeR or DESeq
model fitting, TargetScan prediction, or historical database enrichment has
been executed. Successful helper tests cannot establish reproduction of the
paper's gene lists, networks, figures, statistical significance or conclusions.
