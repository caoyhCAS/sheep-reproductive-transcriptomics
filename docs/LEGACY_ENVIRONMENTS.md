# Legacy statistical environments and analysis boundaries

These scripts are a **method-based reconstruction, not archived original author code**. The paper names edgeR for RSEM counts from M2's ovary/endometrium libraries and DESeq **1.22.0** for sheep/conserved miRNAs. It does not specify edgeR's version, normalization details, dispersion settings, model formula, or low-count filtering choices. The decisions below are explicit reconstruction assumptions.

## Environments are candidates, not verified installations

| File | Purpose | Limit |
|---|---|---|
| `envs/validation.yml` | Current base R plus Python for contract tests | Does not install or run either statistical package |
| `envs/edger-reconstruction.yml` | Portable edgeR environment | Unpinned edgeR is intentional because the study version is unknown; not a historical lock |
| `envs/deseq-1.22.0.yml` | Requested historical DESeq version | Candidate constraint only; historical binaries/dependencies may no longer resolve |

No environment has been solved here. Installing a package is separate from validating its behavior. Archive the resolved packages, `sessionInfo()`, source tarball checksums, platform and compiler details. `edger_denovo.R --expected-version VERSION` can enforce the version you chose. `deseq_mirna.R` requires `packageVersion("DESeq") == "1.22.0"`; it rejects a different patch release, a missing package, and any attempted silent DESeq2 replacement.

DESeq 1.22 belongs to the Bioconductor 3.2 era. The surviving [RELEASE_3_2 branch](https://github.com/bioconductor-source/DESeq/blob/RELEASE_3_2/DESCRIPTION) identifies 1.22.1, so checking out that branch is **not** proof that you installed the paper's 1.22.0. Obtain and verify the exact historical source and compatible dependencies separately. The scripts do not install packages or download executable code. The package APIs were checked against the official-source mirror; the exact historical numerical implementation has not been executed in this reconstruction.

The mapping/assembly versions reported in the paper include FastQC 0.11.4, Cutadapt 1.9.1, TopHat 2.0.8b, Cufflinks/Cuffmerge/Cuffdiff 2.1.1, Trinity 2.1.1, PRINSEQ-LITE 0.20.4, miRDeep2 0.0.7, Blast 2.4.0+, TargetScan 7.0 and Cytoscape 3.4.0. RSEM's version is unreported. Use separately reviewed environments for these tools. Modern aligners, DESeq2 or current annotation databases are not automatic substitutes. The apparent Bowtie/Bowtie2 and reference-assembly discrepancies need explicit resolution before running those stages.

## Input contracts

Both R scripts take `--counts`, `--samples`, `--numerator`, `--denominator`, and `--output`. Inputs are tab-separated, with literal headers:

- Counts: `gene_id` followed by one column per sample. IDs are unique; all counts are finite and nonnegative. `gene_id` is the shared downstream identifier even for miRNAs.
- Samples: `sample_id`, `condition`, `individual_id`. Sample IDs must exactly match count columns; order is aligned by ID. The script selects the two named conditions. Prepare a comparison-specific table from the study sample manifest; do not silently create independent-animal IDs for repeated libraries.
- DESeq additionally requires `--features`: `gene_id`, `class`, with class exactly `sheep`, `conserved`, or `novel`. The feature set must match the count matrix. Only nonzero sheep/conserved rows are fitted; novel and all-zero rows remain `NOTEST` in the export. This normalization/test feature universe is a reconstruction choice, not a recovered author setting.

edgeR accepts fractional RSEM expected counts without rounding. They must be expected fragment counts, not TPM or FPKM. DESeq accepts integer read counts only and rejects fractions instead of rounding them. Its input constructor is documented for integer counts in the [legacy source](https://github.com/bioconductor-source/DESeq/blob/RELEASE_3_2/man/newCountDataSet.Rd). No pseudocount is added to rescue size-factor estimation.

`--validate-only` checks these input/design contracts using base R, creates no output files, and deliberately bypasses package loading. It cannot establish that an analysis environment is installed or numerically correct. Without this flag, an absent/wrong required package fails before data processing.

## Biological independence is a hard boundary

The M2 de novo comparison uses M2-OA, M2-OB, M2-EA and M2-EB: **all four libraries come from one animal**. Mouflon ovarian libraries also include repeated material from each animal; endometrium has only M2 among mouflons. These libraries do not create independent biological replicates. The exact tests implemented here have no individual blocking term and do not estimate a paired or mixed model.

Both scripts reject any comparison containing repeated `individual_id` values or fewer than two distinct individuals per condition unless `--allow-pseudoreplication` is supplied. That flag authorizes an explicitly labeled **descriptive historical calculation**; it does not make the P values valid for population inference. The result's `analysis_scope` and provenance preserve this limitation. A redesigned analysis that aggregates suitable technical replicates or models biological individuals is a separate analysis, not a silent replacement of the historical comparison.

If edgeR has only one library per condition, it additionally requires a caller-justified `--dispersion NUMBER`; no biological dispersion is invented. If DESeq has no replicated condition, only an explicit `--dispersion-method blind` is accepted. This exploratory choice ignores condition labels during dispersion estimation and is not evidence that the paper used it.

## Explicit analysis choices

`edger_denovo.R` uses TMM library normalization, classical common then tagwise dispersion estimation, and `exactTest`. The denominator is the first group and numerator the second. Its double-tail test, `prior.count=0.125`, and `big.count=900` are reconstruction choices from the [edgeR API](https://github.com/bioconductor-source/edgeR/blob/RELEASE_3_4/man/exactTest.Rd); the remaining dispersion defaults come from the installed version and are recorded as such. All-zero features are not tested; no additional low-expression filter is silently introduced. P values are adjusted by Benjamini–Hochberg over the tested features. The paper's de novo threshold, `abs(log2fc) >= 2` and `padj <= 0.01`, is applied downstream.

```bash
Rscript scripts/edger_denovo.R --counts data/M2.rsem_counts.tsv --samples data/M2.contrast.tsv --numerator ovary --denominator endometrium --output results/M2.edger.tsv --allow-pseudoreplication
```

DESeq requires explicit `--dispersion-method pooled|per-condition|blind`, `--sharing-mode maximum|fit-only|gene-est-only`, and `--fit-type parametric|local`. These choices are not reported in the paper; the example selects one conventional option set without claiming historical equivalence. The [legacy dispersion documentation](https://github.com/bioconductor-source/DESeq/blob/RELEASE_3_2/man/estimateDispersions.Rd) defines their different assumptions. The script performs `estimateSizeFactors`, `estimateDispersions`, then `nbinomTest(condA=denominator, condB=numerator)`, whose [documented ratio is B/A](https://github.com/bioconductor-source/DESeq/blob/RELEASE_3_2/man/nbinomTest.Rd). The paper's miRNA threshold is `abs(log2fc) >= 1` and `padj <= 0.05`, applied downstream.

```bash
Rscript scripts/deseq_mirna.R --counts data/mirna_counts.tsv --samples data/mirna.contrast.tsv --features data/mirna_classes.tsv --numerator European_mouflon_ovary --denominator Finnsheep_ovary --output results/mirna.deseq.tsv --dispersion-method pooled --sharing-mode maximum --fit-type parametric --allow-pseudoreplication
```

Each successful analysis writes the requested TSV, `OUTPUT.provenance.tsv`, and `OUTPUT.sessionInfo.txt`; existing files are refused. The standardized columns are `gene_id`, `log2fc`, `padj`, `status`, `numerator`, `denominator`, with `class` added for miRNA. Extra columns include P value, mean expression and analysis scope. Mean expression is normalized CPM for edgeR and normalized base mean for DESeq; these are different units. `status=OK` requires a nonmissing log2 fold change and finite adjusted P value; meaningful positive/negative infinity is retained. Untested/invalid rows have `status=NOTEST` and `padj=NA` and are ignored by downstream selection.

Contract tests use base R to exercise validation, sample matching, count types, contrast orientation arguments, and the independence gate. They do not establish package installation, numerical equivalence, biological independence, or reproduction of the paper's gene lists. When R is unavailable, those tests are explicitly skipped; a successful Python-only test run does not mean the R analyses ran.
