# Study methods and reconstruction scope

This repository organizes the analyses reported by Yang et al. in *Comparative mRNA and miRNA expression in European mouflon (Ovis musimon) and sheep (Ovis aries) provides novel insights into the genetic mechanisms for female reproductive success*, [DOI:10.1038/s41437-018-0090-1](https://doi.org/10.1038/s41437-018-0090-1). It is a method-based reconstruction. The reviewed article and supplementary materials contain workflow diagrams, named upstream programs and results, but no archive of the authors' original analysis scripts or complete command lines. Repository code must therefore not be described as the original study code or as a completed reproduction of its biological results.

The article appeared online on 21 May 2018 and in *Heredity* 122, 172–186 in February 2019. See [SOURCES.md](SOURCES.md) for primary materials and [PARAMETERS.md](PARAMETERS.md) for thresholds and unresolved settings. Executable coverage is described separately in [WORKFLOW.md](WORKFLOW.md).

## Samples and experimental design

The study compared European mouflon with the domestic Finnsheep breed, using ovary and endometrium. mRNA libraries were sequenced as paired-end 100-base reads; miRNA libraries used single-end 50-base reads on Illumina HiSeq 2000.

| Species and tissue | Sample labels | Distinct animals | mRNA libraries | miRNA libraries |
| --- | --- | ---: | ---: | ---: |
| Finnsheep ovary | F1-O through F6-O | 6 | 6 | 6 |
| Finnsheep endometrium | F2-E, F3-E | 2 | 2 | 2 |
| European mouflon ovary | M1-OA, M1-OB, M2-OA, M2-OB, M3-OA, M3-OB | 3 | 6 | 6 |
| European mouflon endometrium | M2-EA, M2-EB | 1 | 2 | 1, M2-EB only |

The 16 mRNA and 15 miRNA libraries are not 31 independent animals. A/B mouflon labels are samples from the same animal; M2-EA and M2-EB were collected from the two uterine horns of one ewe. Finnsheep ovarian samples were collected during follicular growth and endometrial samples during early pregnancy. Mouflon cycle stages were not known, except that M2 had a corpus luteum. Preserve animal, tissue and library identifiers separately. Species, reproductive stage, sampling conditions and within-animal sampling can affect expression comparisons.

For inference at the animal or species level, treating A/B samples as unrelated biological replicates is inappropriate. In particular, the mouflon endometrium has only one animal, and only one miRNA library. A reanalysis cannot manufacture biological replication. The historical pairwise contrasts can be documented or run with explicitly supplied statistical assumptions, but their sampling limitations remain.

## Reference-based mRNA analysis

1. Inspect raw reads with FastQC 0.11.4 and remove adapters with Cutadapt 1.9.1.
2. Align with TopHat 2.0.8b. The article names Oar v4.0 and Ensembl release 83 annotation together; these need resolution before execution because the official release 83 sheep GTF is for Oar v3.1.
3. Assemble transcripts with Cufflinks 2.1.1, then use Cuffmerge to construct merged annotation sets. The paper describes merging ovarian or endometrial samples within each species, but does not archive the final GTF supplied to every cross-species comparison.
4. Use Cuffdiff from Cufflinks 2.1.1 for expression and differential-expression analysis. The study reports FPKM and selects genes with absolute log2 fold change at least 2 and adjusted P value at most 0.01.
5. Submit candidate genes to DAVID 6.8 for GO/pathway enrichment, using adjusted P value below 0.05. The background gene universe, identifier conversion and exact correction method for this DAVID branch are not specified. The paper created an enrichment word cloud with the [Jason Davies online word-cloud generator](https://www.jasondavies.com/wordcloud/), without supplying plotting settings.

Comparisons include mouflon versus Finnsheep within ovary, mouflon versus Finnsheep within endometrium, and ovary versus endometrium within mouflon. A valid execution must provide the exact sample grouping, shared coordinate system and comparison orientation rather than infer these from file order. Cuffdiff output should be filtered on its reported q value and successful test status; expression values alone do not supply a differential-expression test.

## De novo mouflon transcriptome

All eight mouflon mRNA libraries were assembled with Trinity 2.1.1. Assembly assessment included read remapping. RSEM estimated transcript abundances for the four M2 libraries, and edgeR supplied differential-expression tests. Reported thresholds are absolute log2 fold change at least 2 and FDR at most 0.01.

Supplementary Figure S2 names Trinity's assembly and supporting abundance/matrix/differential-expression scripts. Its labels include `Trinity.pl`, `Align_and_estimate_abundance.pl`, `Abundance_estimate_to_matrix.pl`, and `Run_DE_analysis.pl`; labels in the diagram are not guaranteed to be the case-sensitive filenames of an installed release. These are upstream tool components, not an archived collection of author-written scripts.

The Results list four individual-library contrasts: M2-OA versus M2-EA, M2-OA versus M2-EB, M2-OB versus M2-EA and M2-OB versus M2-EB. edgeR dispersion settings for these single-library contrasts, RSEM/edgeR versions, transcript-to-unigene reduction and normalization options are not reported. They must be supplied and recorded rather than hidden in a default.

Differentially expressed sequences were searched against NCBI nt with BLASTN 2.4.0+, an E-value threshold of 1e-10 and a reported 100% matching criterion. The article does not define whether that criterion means identity, query coverage or both. PANTHER 10 provided enrichment, with Bonferroni-adjusted P value below 0.05. Preserve the exact nt/database snapshot and interpretation of the matching criterion.

## miRNA analysis

FastQC and Cutadapt were followed by PRINSEQ-LITE 0.20.4 for further processing, including ambiguous-read filtering, reformatting and adapter trimming. The ambiguous-base threshold and exact division of trimming between the programs are unspecified. Reads of 18–26 nucleotides were retained. The Methods specify an Oar v3.1 index and a program called Bowtie version 2.2.1; the name, version and cited Bowtie paper are inconsistent. Supplementary Figure S1 explicitly includes the miRDeep2 mapper and prediction scripts. Resolve the installed aligner and mapping format before connecting this branch to miRDeep2 0.0.7.

Predicted miRNAs were compared with miRBase 21 and grouped as sheep-known, conserved in other mammals, or novel. Conserved matches allow fewer than four mismatches. The prose combines a reads-at-most-10 condition and a score-below-5 condition without an unambiguous Boolean filtering rule. See [PARAMETERS.md](PARAMETERS.md).

DESeq 1.22.0, the historical DESeq package rather than DESeq2, was used for differential expression of sheep-known and conserved miRNAs. The stated thresholds are absolute log2 fold change at least 1 and adjusted P value at most 0.05. Novel predictions were characterized separately and are not described as inputs to that differential-expression/network branch. Count normalization, dispersion fitting and treatment of the unequal biological replication are not specified. Use original integer counts for the model, not FPKM values from the mRNA branch.

## Integrated miRNA–mRNA network

TargetScan 7.0 predicted miRNA targets from seed complementarity to 7-mer and 8-mer sites in mRNA 3-prime UTRs. Predicted targets were intersected with mRNA DEGs, then restricted to miRNA–mRNA pairs with opposing expression directions and context++ score percentiles of at least 50. The article's phrase about inverse correlation does not provide a Pearson/Spearman coefficient, significance threshold, matched-sample procedure or alternative statistical test; a direction-based intersection should not be presented as a measured correlation.

All full networks were visualized with Cytoscape 3.4.0. A main network retained the top ten downregulated miRNAs. The main figure identifies Finnsheep-downregulated/mouflon-upregulated miRNAs, while other networks reverse the comparison. The rank variable, tie policy and whether selection precedes target filtering are not completely specified. Preserve the comparison orientation and expose the ranking rule. The TargetScan UTR sequence collection, transcript/species mapping, conservation settings, isoform aggregation and context++ prediction files are necessary inputs that the paper does not archive.

DAVID 6.8 supplied target-gene GO/pathway enrichment. miRTarBase and prior literature were checked for support of predicted interactions; this is distinct from experimentally validating every edge. The study also compared mRNA candidates with the Sheep QTL database and drew heatmaps using pheatmap 1.0.8. Historical database versions, plotting normalization, clustering distance and complete plotting scripts are unavailable.

## Practical reproduction limits

The public read project provides 31 runs, and its experiment aliases identify the paper's libraries. The project's shared sample alias is a list reused across all runs, so it is unsuitable as an individual-library identifier. Raw reads are approximately 212 billion bases in the current archive; repository smoke tests and synthetic fixtures do not analyze this data volume or recreate the paper's results.

Before scientific execution, resolve the reference/annotation pair, alignment version, adapters and library strandedness, biological grouping, differential-expression model assumptions, miRNA precursor filter, TargetScan inputs and enrichment background. Tool substitutions and new statistical decisions should be recorded as a reanalysis, with differences from the historical method made explicit.
