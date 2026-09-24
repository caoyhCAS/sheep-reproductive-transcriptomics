# Reported parameters and unresolved choices

The values below come from the [article](https://doi.org/10.1038/s41437-018-0090-1), its Methods, main figure captions and Supplementary Figures S1–S2. They describe the study rather than guarantee installed versions or complete executable defaults. Tool invocation details are in [WORKFLOW.md](WORKFLOW.md).

## Reported settings

| Analysis | Software or resource | Setting reported |
| --- | --- | --- |
| Raw-read QC | FastQC 0.11.4 | Both RNA types; no filtering thresholds supplied |
| Adapter removal | Cutadapt 1.9.1 | Both RNA types; adapter sequences/options unspecified |
| mRNA alignment | TopHat 2.0.8b | Oar v4.0 with Ensembl 83 annotation claimed; requires resolution |
| mRNA quantification | Cufflinks/Cuffmerge/Cuffdiff 2.1.1 | FPKM; `abs(log2FC) >= 2`, adjusted P `<= 0.01` |
| De novo assembly | Trinity 2.1.1 | All eight mouflon mRNA libraries |
| De novo quantification | RSEM, version unspecified | M2-OA, M2-OB, M2-EA, M2-EB |
| De novo expression | edgeR, version unspecified | `abs(log2FC) >= 2`, FDR `<= 0.01` |
| De novo annotation | BLASTN 2.4.0+, NCBI nt | E value `1e-10`; reported matching ratio 100% |
| De novo enrichment | PANTHER 10 | Bonferroni-adjusted P `< 0.05` |
| miRNA preprocessing | PRINSEQ-LITE 0.20.4 | Retained length `18 <= nt <= 26` |
| miRNA alignment | Bowtie 2.2.1 as written | Oar v3.1/Ensembl 83; program/version conflict |
| miRNA prediction | miRDeep2 0.0.7 | Removal criteria mention reads `<= 10` and score `< 5` |
| miRNA classification | miRBase 21 | Sheep-known, conserved or novel; conserved mismatch count `< 4` |
| miRNA expression | DESeq 1.22.0 | Known/conserved miRNAs; `abs(log2FC) >= 1`, adjusted P `<= 0.05` |
| Target prediction | TargetScan 7.0 | Seed positions 2–8; 7-/8-mer UTR sites; context++ score percentile `>= 50` |
| Network integration | Expression/target intersection | Opposite miRNA/mRNA expression directions |
| Network display | Cytoscape 3.4.0 | Full networks and a top-ten-downregulated-miRNA main network |
| GO/pathway enrichment | DAVID 6.8 | Adjusted P `< 0.05` |
| Heatmaps | pheatmap 1.0.8 | Mentioned in Figure 2 caption; complete plotting options unspecified |

The main text uses inclusive differential-expression cutoffs. Table 2's miRNA footnotes use strict inequalities. A reconstruction should choose and document its boundary convention; the main-text inclusive thresholds above are not a claim that both descriptions are identical.

## Choices that must not be silently filled in

### Assembly and annotation

The mRNA Methods name Oar v4.0 and Ensembl release 83. The [official release 83 GTF directory](https://ftp.ensembl.org/pub/release-83/gtf/ovis_aries/) identifies `Ovis_aries.Oar_v3.1.83.gtf.gz`. The miRNA Methods explicitly name Oar v3.1, but the Results describe both RNA types as mapped to Oar v4.0. No coordinate conversion or replacement annotation is given. A v3.1 GTF must not simply be applied to a v4.0 FASTA. Use a confirmed compatible pair, retain checksums/assembly identity and label departures from the paper.

### Bowtie and miRDeep2

The paper writes Bowtie 2.2.1 but cites the original Bowtie publication and says Bowtie-build. Version 2.2.1 belongs to Bowtie 2; the two programs have distinct executable names and mapping conventions. Index and output compatibility must be checked for the selected historical versions, along with the mapper/prediction handoff for the chosen miRDeep2 release. Do not turn an unresolved name into an apparently pinned historical command.

### Precursor filters

Two interpretations of the reads/score prose are possible:

- Joint rejection: discard only candidates with `reads <= 10 AND score < 5`.
- Independent quality requirements: retain only candidates with `reads > 10 AND score >= 5`, which discards when either requirement fails.

The difference affects candidates with high read counts but low scores, or low counts but high scores. Require a named filter mode or an explicit recorded reconstruction decision. The equality boundaries are also meaningful: ten reads are in the low-read category, while a score of five is not below the score cutoff.

### Count models and biological replication

The paper does not supply DESeq/edgeR normalization, dispersion, design matrix, replicate-handling, low-count filtering or seed settings. Within-animal A/B libraries cannot be treated as independent animals. The de novo Results enumerate four single-library contrasts without reporting the dispersion used for their edgeR tests. A fabricated dispersion or a replacement DESeq2 model would be a new analysis. Known/conserved-versus-novel inclusion and the direction of each comparison must be explicit.

### BLAST matching and unigene annotation

The 100% matching ratio is not defined as percent identity, alignment coverage or both. Preserve `pident`, aligned length, query length and query coverage in tabular outputs so a user can apply a documented interpretation. E-value filtering alone does not implement this criterion. The nt snapshot, best-hit/tie rule, transcript-to-gene mapping and unigene representative rule remain unspecified.

### TargetScan and network construction

TargetScan requires compatible miRNA family/seed definitions and UTR sequences. The paper does not supply the sheep UTR collection, genome build, transcript identifiers, other-species/conservation alignments, isoform selection rule or context++ tables. A percentile is not the raw context++ score: do not filter raw scores at 50. Preserve transcript-to-gene mappings and disclose whether multiple sites/isoforms are retained, summed or reduced to a best score.

“Top ten downregulated” lacks a fully specified sorting variable and tie policy. Do not silently sort by adjusted P value, absolute fold change or target degree as if the paper required it. Select within a stated contrast and tissue, record the direction, define ties, and state whether ranking happens before or after target intersection. Positive numbers in a table organized into separate up/down panels do not establish a universal fold-change sign.

### Other unreported inputs

Adapter sequences, quality trimming thresholds, library strandedness, TopHat intron/mismatch settings, multimapping treatment, merge strategy across species, mature-miRNA duplicate handling, miRBase species mapping, database release dates, DAVID/PANTHER backgrounds, GO correction details outside the explicit PANTHER branch, and heatmap normalization/clustering choices are not fully specified. They should remain documented inputs or clearly identified implementation choices.
