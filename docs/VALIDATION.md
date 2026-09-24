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

Remote upload and GitHub CI have **not yet run** for this source package.
The target repository must exist before the connected uploader can publish it.
After publication, the commit's actual Actions result is authoritative; this
file records the pre-publication local validation state.

No TopHat/Cufflinks/Trinity/RSEM/miRDeep2 biological computation, edgeR or DESeq
model fitting, TargetScan prediction, or historical database enrichment has
been executed. Successful helper tests cannot establish reproduction of the
paper's gene lists, networks, figures, statistical significance or conclusions.
