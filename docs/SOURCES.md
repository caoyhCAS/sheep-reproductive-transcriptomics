# Source provenance

Primary materials were inspected on **2026-09-24**. This repository is a method-based reconstruction rather than an author-supplied historical code archive. The primary paper and its supplementary materials describe the workflow, but the reviewed files contain no complete original script package. Targeted searches by DOI, title and SRA project did not locate such an archive; this does not establish that private or subsequently deposited author code does not exist.

## Article

Yang J, Li X, Cao Y-H, Pokharel K, Hu X-J, Chen Z-H, Xu S-S, Peippo J, Honkatukia M, Kantanen J, Li M-H. Comparative mRNA and miRNA expression in European mouflon (*Ovis musimon*) and sheep (*Ovis aries*) provides novel insights into the genetic mechanisms for female reproductive success. *Heredity* **122**, 172–186 (2019). Published online 21 May 2018.

- DOI and publisher: <https://doi.org/10.1038/s41437-018-0090-1>
- Publisher article: <https://www.nature.com/articles/s41437-018-0090-1>
- PubMed record: <https://pubmed.ncbi.nlm.nih.gov/29784930/>
- PMC record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6327046/>
- Authors' institutional repository record: <https://jukuri.luke.fi/items/755c07a5-ca8c-404c-9d6c-4f4d2ccb5c6e>
- Institutional publisher-PDF mirror used for full-text extraction: <https://jukuri.luke.fi/server/api/core/bitstreams/18c58fcb-8391-45d1-a19c-e0fedd294515/content>

The article is licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), subject to its stated exception for separately credited third-party material. New repository implementation code has its own license; the article license does not automatically relicense external software.

Publisher redirects and a failure of the Europe PMC XML endpoint prevented direct XML retrieval during the audit. The institutional publisher PDF, public repository metadata and publisher-indexed text provided the article evidence. The Europe PMC supplementary-files endpoint succeeded.

## Supplementary material

Public archive endpoint: <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6327046/supplementaryFiles>

| Original filename | Description | SHA-256 |
| --- | --- | --- |
| `41437_2018_90_MOESM1_ESM.xlsx` | Supplementary Table S11, endometrial mRNA differential-expression results, per archive label | `681c8329f2d36798c134b71ea8f75c12481fc35be6175a860e1406cef6f35365` |
| `41437_2018_90_MOESM2_ESM.doc` | Supplementary Results, Figures S1–S5, Tables S1–S23 and references; S10/S11 supplied separately | `42c2c216a90028226f39d79a2a2c85e6173361eb1b5c685d67f0c4bc7952c233` |
| `41437_2018_90_MOESM3_ESM.xlsx` | Supplementary Table S10, ovarian mRNA differential-expression results, per archive label | `53a3eee9188a6408295d9fe1eeebb0ba964855a639b9ed9e6cc53ea6dd49045d` |

The primary PDF used for extraction has SHA-256 `0748b1057efdb425c4aebb8294695d9f5c40a9ab60479a13259bb48b37d9a31a`.

The spreadsheet descriptions above follow the public archive labels; their cell-level differential-expression results were not independently reanalyzed during the methods audit.

The DOC was converted locally for read-only inspection. Its Figure S1 identifies the reference-mRNA, miRNA and integration paths. Figure S2 identifies the Trinity de novo path. The figures were visually checked because workflow labels are embedded in images and absent from plain-text extraction. These diagrams name upstream scripts, but do not provide their command lines or custom source code.

## Public sequencing data

- BioProject: [PRJNA451237](https://www.ncbi.nlm.nih.gov/bioproject/451237)
- SRA study: [SRP142554](https://www.ncbi.nlm.nih.gov/sra/?term=SRP142554)
- ENA project: <https://www.ebi.ac.uk/ena/browser/view/PRJNA451237>
- ENA machine-readable run metadata: <https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA451237&result=read_run&format=json>

The project has 31 experiments/runs: 16 paired RNA-Seq and 15 single-end miRNA-Seq libraries. The experiment aliases/library names identify the sample labels. The shared BioSample/sample-alias field repeats a multi-sample list, so it should not be used as the sample mapping. Download URLs and MD5 checksums should be retrieved from ENA run metadata and verified after download. Raw data and external references are not included in this repository.

## References and external programs

The study's precise versions and uses are cataloged in [PARAMETERS.md](PARAMETERS.md). External tools and databases remain subject to their own terms and licenses. Install them from their maintainers rather than assume they are covered by the repository license.

| Resource | Primary source or project |
| --- | --- |
| Ensembl release 83 sheep GTF | <https://ftp.ensembl.org/pub/release-83/gtf/ovis_aries/> |
| Ensembl release 83 sheep genomic FASTA | <https://ftp.ensembl.org/pub/release-83/fasta/ovis_aries/dna/> |
| Trinity and its bundled utilities | <https://github.com/trinityrnaseq/trinityrnaseq> |
| miRDeep2 and its mapper/prediction utilities | <https://github.com/rajewsky-lab/mirdeep2> |
| Bowtie 2 release history | <https://bowtie-bio.sourceforge.net/bowtie2/news.shtml> |
| miRBase | <https://www.mirbase.org/>; the paper specifies release 21 |
| TargetScan | <https://www.targetscan.org/>; the paper specifies version 7.0 |
| DESeq | <https://bioconductor.org/packages/DESeq/>; historical version 1.22.0 |
| edgeR | <https://bioconductor.org/packages/edgeR/>; historical version not specified |
| DAVID | <https://david.ncifcrf.gov/>; historical version 6.8 |
| PANTHER | <https://www.pantherdb.org/>; historical version 10 |
| Cytoscape | <https://cytoscape.org/>; historical version 3.4.0 |
| Sheep QTL database | <https://www.animalgenome.org/cgi-bin/QTLdb/OA/> |

Current project landing pages are discovery links, not evidence that a current release or database snapshot reproduces a historical result. The exact author-used releases, UTR/annotation files and database snapshots should be preserved with checksums when available. No third-party script is claimed to be original author code merely because it appears in a supplementary workflow diagram.
