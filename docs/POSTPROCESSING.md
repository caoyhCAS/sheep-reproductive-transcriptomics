# Export contracts and downstream filters

These scripts are newly written method reconstructions, not original author
scripts. Run each command with `--help`. They consume explicit exports; none
substitutes for an external statistical model or database search.

## Differential expression and orientation

`scripts/expression_analysis.py normalize-cuffdiff` reads `gene_exp.diff` with
`gene_id`, `sample_1`, `sample_2`, `status`, `value_1`, `value_2`,
`log2(fold_change)`, and `q_value`. It selects the requested denominator as
sample_1 and numerator as sample_2; it preserves the native fold change to
avoid changing boundary decisions through rounded FPKM values. Only `OK`
results enter filtering. One-zero infinite fold changes can be retained;
NaN fold changes and invalid adjusted P values in tested rows are errors.

Normalized exports require these tab-separated columns:

```text
gene_id  log2fc  padj  status  numerator  denominator
```

Here `gene_id` means the unique feature ID, including a mature miRNA ID in a
miRNA table. miRNA tables additionally require `class`: `sheep`, `conserved`,
or `novel`. Only sheep and conserved miRNAs enter differential-expression
filtering. The R scripts generate this contract. `filter --mode mrna` and
`denovo` require absolute log2FC >=2 and adjusted P <=0.01; `mirna` requires
absolute log2FC >=1 and adjusted P <=0.05.

An additional `analysis_scope` column is preserved from R outputs through
filtering; it defaults to `unspecified` if absent. Network files carry separate
`mrna_analysis_scope` and `mirna_analysis_scope` columns, preserving the warning
`descriptive_nonindependent_libraries` where relevant. For historical Cuffdiff
comparisons, supply `--analysis-scope descriptive_nonindependent_libraries`
when normalizing; the normalizer cannot infer animal identities from gene_exp.diff.

## Preparing RSEM expected-count matrices

`scripts/assemble_counts.py` joins per-sample RSEM exports by feature ID and
retains fractional `expected_count` values. Its manifest has `sample_id` and
`rsem_results` columns, with result paths relative to the manifest. Select
`--level genes` for `.genes.results` or `--level isoforms` for
`.isoforms.results` explicitly; both export the generic feature key `gene_id`.
The term “unigenes” in the paper does not recover the authors' exact aggregation
choice. The helper rejects mismatched feature sets, duplicate IDs, missing
values, reused input files, and attempts to overwrite existing outputs.

```bash
python scripts/assemble_counts.py --manifest data/M2.rsem_manifest.tsv --level genes --output results/M2.counts.tsv
```

For the reported M2 de novo comparison include M2-OA, M2-OB, M2-EA, M2-EB and
record that all four libraries belong to M2 in the R comparison sample table.
Do not feed TPM or FPKM to a count model.

## TargetScan and Cytoscape

Provide an externally computed TargetScan 7.0 export normalized to:

```text
mirna  gene_id  context_percentile
```

Percentiles are 0..100, not the raw context++ score. IDs must match the
expression files exactly. Resolve gene/transcript/family mappings before
export, retaining a mapping audit. The script does not infer missing UTR
alignments, conservation settings, or context++ model inputs.

`integrate` requires matching numerator/denominator across mRNA and miRNA,
reapplies both DE thresholds, retains percentiles >=50, and requires opposite
nonzero fold-change signs. Multiple exported rows for a pair are collapsed by
maximum percentile (a reconstruction decision). Network output is a TSV edge
list for Cytoscape import, not a reproduction of the original layout. The
core file selects the 10 most downregulated network-participating miRNAs by
log2FC, breaking ties by ID; the paper leaves the precise ranking convention
unspecified. Opposite expression directions and target predictions do not
establish causality or experimental validation.

## Precursor, annotation, and enrichment filters

`scripts/filter_evidence.py` provides three explicit export filters:

| Command | Required TSV columns | Decision |
|---|---|---|
| `precursors` | precursor_id, reads, score | `--policy joint-failure` rejects reads <=10 **and** score <5; `either-failure` rejects either criterion. Caller must choose. |
| `blast` | qseqid, sseqid, evalue, pident, qcovhsp | E-value <=1e-10 plus `--matching-ratio identity`, `coverage`, or `both` at 100%. Coverage is per HSP. |
| `enrichment` | term_id, adjusted_p, correction | Adjusted P <0.05, and correction must match the explicit `--correction` label. |

Precursor read/score exports must first be mapped explicitly from miRDeep2's
result columns. No mature-miRNA sequence matcher is silently invented: classify
against miRBase21 externally, document alignment/gap/mismatch conventions,
and preserve sheep/conserved/novel assignments. The reported conserved criterion
is fewer than four mismatches, but alignment details are unreported.

BLAST exports can be created using BLAST+2.4.0 with `-outfmt '6 qseqid sseqid
evalue pident qcovhsp'`, adding that TSV header before filtering. The original
meaning of “matching ratio 100%” is unresolved; choices have different results.
Provide the NCBI nt snapshot and checksum rather than treating today's database
as the historical one. This helper does not select one best ortholog from
multiple hits.

PANTHER10 used Bonferroni correction. For DAVID6.8, the exact adjusted-P column
and background universe were not specified: record both when exporting, do
not relabel unadjusted P as adjusted P. The helper only filters supplied
corrected values; it does not run enrichment. Historical web services,
miRTarBase evidence lookup, and Wordcloud/Cytoscape layout remain external
steps with explicitly documented inputs and version limitations.
