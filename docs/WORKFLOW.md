# Legacy workflow reconstruction

This repository contains newly written reconstruction code for
`10.1038/s41437-018-0090-1`. It is not the authors' original pipeline. A plan can
be reviewed without installing any biological tools or obtaining sequence data.
Execution uses the named legacy tools; it does not replace TopHat/Cufflinks,
DESeq or other historical methods with modern alternatives.

## Plan first

From the repository root:

```bash
python scripts/workflow.py --config config/paper.json --samples config/samples.tsv --stage all --output results/plan.json
```

This creates only the requested JSON plan. Without `--output`, JSON goes to
standard output. Missing adapters, unresolved scientific choices and unreviewed
acknowledgments appear in `execution_blockers`; they do not prevent inspecting
the plan. `<UNRESOLVED:...>` strings are deliberately nonexecutable placeholders.
Planning never downloads sequence data or invokes a tool.

| Stage | Included steps | Main products |
| --- | --- | --- |
| `mrna` | FastQC, paired Cutadapt, TopHat, Cufflinks for all 16 libraries | Trimmed FASTQs, `accepted_hits.bam`, per-library `transcripts.gtf` |
| `cuffdiff` | `mrna`, four species/tissue Cuffmerge groups, resolved transcriptome selection, three Cuffdiff contrasts | `merged.gtf`, `gene_exp.diff` |
| `trinity` | FastQC/Cutadapt for all eight mouflon libraries; Trinity assembly | `Trinity.fasta` |
| `rsem` | Mouflon read preprocessing, an existing Trinity FASTA, RSEM gene map/reference and four M2 quantifications | `.genes.results`, `.isoforms.results` |
| `denovo` | `trinity` followed by the RSEM steps | Assembly and M2 quantifications |
| `mirna` | FastQC, single-end Cutadapt, PRINSEQ for 15 libraries | Clean 18–26 nt FASTQ files |
| `all` | All automated stages above | Their combined products |

All paths are relative to the working directory. Execution deliberately accepts
only ASCII letters/digits, underscores, dots, slashes and hyphens in data/output
paths and the working directory. Legacy Cuffmerge and Trinity construct shell
commands internally, so argument lists at the Python boundary alone would not
make arbitrary filenames safe. Stage files under safe paths before using these
archived programs. Paths cannot begin with `-`. FASTQ paths must be unique, and
paired mate filenames must have distinct basenames for FastQC output names.

## Sample and biological-individual contract

`config/samples.tsv` has exactly the 16 reported sample IDs. Its columns are
`sample_id`, `species`, `tissue`, `individual_id`, `tissue_replicate`, `mrna_r1`,
`mrna_r2`, and `mirna_fastq`. The FASTQ paths are placeholders; replace them
with downloaded/verified files. These metadata are checked against the reported
design rather than inferred from a shared archive BioSample accession.

- Finnsheep ovary: F1-O through F6-O, six different ewes.
- Finnsheep endometrium: F2-E and F3-E, from F2 and F3.
- European mouflon ovary: M1-OA/OB, M2-OA/OB and M3-OA/OB, six tissue samples
  from **three animals**.
- European mouflon endometrium: M2-EA and M2-EB, both from **one animal, M2**.
- M2-EA has no miRNA library. Its `mirna_fastq` must remain `.`; the other 15
  samples each have one single-end miRNA FASTQ.

The A/B labels are tissue replicates from the same animal, not new biological
individuals. All four RSEM quantifications are also from M2. Reconstructed
historical contrasts can therefore involve pseudoreplication, confounded tissue
sampling and an endometrial group with only one wild individual. Acknowledging
this limitation does not make the design biologically replicated. Do not treat
these reconstructed contrasts as independent evidence of population-level
species differences.

## Choices required before execution

Copy `config/paper.json` to a local configuration and resolve its `null` values
for the selected stage. Computational resources `threads=4` and
`trinity.max_memory=32G` are editable execution examples, not reported paper
parameters. Actual data at the paper's scale can need substantially more memory.

| Field | Required decision |
| --- | --- |
| `mrna.adapter_r1`, `mrna.adapter_r2`, `mrna.adapter_resolution_note` | Actual IUPAC adapter sequences and their provenance; no TruSeq sequence is assumed |
| `mrna.library_type` | `fr-unstranded`, `fr-firststrand` or `fr-secondstrand`; not reported precisely |
| `mrna.cufflinks_annotation_mode` | `guided` (`-g`) or `unguided`; the assembly-guidance option is not specified |
| `reference.genome_fasta`, `annotation_gtf`, `genome_assembly`, `annotation_assembly`, `resolution_note` | A verified mutually compatible FASTA/GTF pair and written resolution of Oar4/Ensembl83 inconsistency |
| `reference.bowtie2_index_prefix`, `bowtie2_index_files` | All six matching `.bt2` or `.bt2l` files for the provided reference; no implicit index build |
| `cuffdiff.transcriptome_policy`, `resolution_note` | Explicit `union_per_contrast` or `provided_per_contrast` interpretation of the merged-GTF ambiguity |
| `cuffdiff.contrast_gtfs` | Three supplied GTF paths when using `provided_per_contrast` |
| `trinity.strand` | `unstranded`, `RF` or `FR`; no hidden strand inference |
| `rsem.aligner`, `strandedness`, `resolution_note` | Explicit `bowtie`/`bowtie2` and `none`/`forward`/`reverse`, because these RSEM settings are unreported |
| `tools.rsem_extract`, `tools.rsem_prepare`, `tools.rsem_quant` | Installed command and actual expected version; the paper does not identify an RSEM version |
| `mirna.adapter`, `adapter_resolution_note` | Actual miRNA adapter sequence/provenance |

The paper describes mapping mRNA to Oar v4.0 with Ensembl release83 annotation;
release83 is associated with Oar v3.1. The workflow refuses a configuration whose
declared FASTA and GTF assembly identifiers differ, even after acknowledgment.
Supply a checked assembly-consistent pair and document any conversion or chosen
departure. Matching declarations alone cannot verify assembly provenance; retain
download URLs, checksums and annotation preparation records. TopHat2 uses its
Bowtie2 backend in this reconstruction; that dependency/version is not fully
specified for mRNA by the paper and must be recorded with the environment.

Cutadapt removes supplied 3′ adapters with `-a` and paired `-A`; no additional
quality/length trimming thresholds are invented. Residual behavior is the
archived tool's default. PRINSEQ uses `-min_len 18 -max_len 26 -ns_max_n 0`;
`-noniupac` also excludes non-ACGTN symbols. This explicit interpretation of
ambiguous-read removal is recorded in the plan. miRNA reads are not deduplicated
because that would discard abundance information.

## Cuffmerge and contrast direction

The paper describes species/tissue-specific merged GTFs but does not identify
which common transcript set Cuffdiff uses across species. Four Cuffmerge calls
reconstruct those groups. `union_per_contrast` adds a second merge of both
group GTFs for each contrast; this is a documented reconstruction decision.
`provided_per_contrast` instead uses the user's three verified GTFs. Neither
choice claims to recover the original unreported transcript universe.

| Contrast key | Numerator (`sample_2`) | Denominator (`sample_1`) |
| --- | --- | --- |
| `ovary_species` | European_mouflon_ovary | Finnsheep_ovary |
| `endometrium_species` | European_mouflon_endometrium | Finnsheep_endometrium |
| `mouflon_tissues` | European_mouflon_ovary | European_mouflon_endometrium |

Cuffdiff receives the denominator BAM group first and numerator group second,
so its `log2(fold_change)` represents `log2(value_2/value_1)`. The JSON contains
these labels, member biological-individual IDs and the direction. Use matching
numerator/denominator labels in downstream exports; never reverse signs based
only on a filename. The mouflon tissue contrast uses the six ovarian and two
endometrial samples in the reference-mapping branch; the de novo/RSEM branch
uses only M2-OA/OB and M2-EA/EB as described for that analysis.

## Trinity and RSEM limits

Trinity2.1.1 receives eight paired mouflon libraries. Its default minimum contig
length/other assembly settings remain defaults; no normalization or extra
Trimmomatic step is added. `--no_version_check` disables its execution-time update
check. The supplementary diagram names Trinity utilities for abundance and DE
analysis. This reconstruction uses the independently documented RSEM interfaces
directly; it does **not** claim to recover the unreported wrapper options.

RSEM's `extract-transcript-to-gene-map-from-trinity` generates its own required
Trinity identifier map, followed by `rsem-prepare-reference` and paired-end
quantification. Mapping Trinity identifiers into gene groups does not reproduce
an undocumented longest-transcript or separate unigene-selection procedure.
Treat the output as Trinity-gene/transcript abundance estimates. Downstream edgeR
and any rounding/normalization decisions are separate, explicitly documented
steps. An `rsem`-only run should use a new output directory and set
`rsem.trinity_fasta` to an existing assembled FASTA.

## miRNA mapping remains a handoff

The Methods name Oar v3.1, `Bowtie-build`, and “Bowtie v2.2.1”, whereas Results
refer to Oar v4.0; the supplement names miRDeep2 `mapper.pl` and `miRDeep2.pl`.
Bowtie and Bowtie2 are different interfaces, and miRDeep2's mapper input/mapping
format cannot be inferred from this description. The automated `mirna` stage
therefore ends at filtered FASTQ. A reviewed mapper/reference choice and valid
miRDeep2 ARF/collapsed-read inputs are still needed. This is an explicit incomplete
handoff, not a simulated mapping result.

The [official miRDeep2 interface documentation](https://github.com/rajewsky-lab/mirdeep2/blob/master/documentation.html)
supports the following **reviewable templates**, using already filtered FASTQ
and actual miRBase21 FASTAs. Run each sample in its own new working directory;
miRDeep2 creates timestamped files there. Replace `/ABS/...` paths. These are
upstream interface templates, not recovered v0.0.7 commands.

```bash
mapper.pl /ABS/sample.clean.fastq -e -h -i -j -m -p /ABS/index/genome -s sample.collapsed.fa -t sample.genome.arf -v
miRDeep2.pl sample.collapsed.fa /ABS/genome.fa sample.genome.arf /ABS/oar.mature.fa /ABS/related.mature.fa /ABS/oar.hairpin.fa
quantifier.pl -p /ABS/oar.hairpin.fa -m /ABS/oar.mature.fa -r sample.collapsed.fa -t oar -y sample
```

Before use, verify v0.0.7 options and dependencies locally, choose and document
the Oar assembly, and explicitly resolve using Bowtie1 `bowtie-build`/`.ebwt`
indexes instead of the paper's inconsistent Bowtie2.2.1 wording. Preserve mapper
abundance suffixes (`_xN`); collapsing is not abundance-discarding deduplication.
No second adapter trim is applied here. Record mapping/quantifier defaults and
the related-species FASTA selection. Export counts with stable feature IDs
across all samples before DESeq; per-sample novel discoveries are not directly
a shared count matrix. Species/conserved classification and precursor/score
filtering remain explicit downstream steps.

## Execute after resolving the plan

```bash
python scripts/workflow.py --config config/local.json --samples config/samples.tsv --stage mrna --execute
```

Execution requires the stage's acknowledgments to be `true`, all external input
files to exist and be nonempty, and version checks to pass. It refuses existing
nonempty output directories. It is a sequential runner with dependency/output
checks, not a resumable scheduler. For reruns choose a new output directory;
never delete prior results just to satisfy a check. No raw sequence files are
downloaded automatically. Plans are inspectable argument arrays, and subprocess
calls do not invoke a shell.

Plan output paths cannot overwrite the input configuration, sample manifest or
an existing plan. With `--execute`, write an optional `--output` plan outside
the run directory; execution itself writes `execution_plan.json` after checks.

Exact reported versions are enforced: FastQC0.11.4, Cutadapt1.9.1,
TopHat2.0.8b, Cufflinks2.1.1, Trinity2.1.1 and PRINSEQ0.20.4. Cuffmerge's own
script reports `merge_cuff_asms v1.0.0`; the preflight records that separately
from the Cufflinks2.1.1 bundle. Trinity2.1.1 `--version` exits with status1 and
prints a `v` prefix; this documented behavior is accepted for that tool only.
Its version-report code contacts update servers even when the normal-run update
check is disabled; planning invokes no version checks. Other tools normally
require version-command exit0. Configured version commands can be adapted to
local packaging while preserving exact reported versions. Related executables
must come from the same verified distribution; a version string is not a binary
checksum or complete environment lock.

Logs, the executed plan and verified version outputs are saved within the run
directory. Compatibility of ancient interpreters, Java, Bowtie dependencies and
external reference content remains the user's installation responsibility.
Tests cover plan structure, sample identities, scientific-choice gates, path
safety and version-check behavior; they do not validate a real legacy tool run
or reproduce the paper's biological results.

Primary command references:

- [Cutadapt1.9 paired trimming](https://cutadapt.readthedocs.io/en/v1.9/guide.html#trimming-paired-end-reads)
- [TopHat manual](https://ccb.jhu.edu/software/tophat/manual.shtml)
- [Cuffmerge](https://cole-trapnell-lab.github.io/cufflinks/cuffmerge/) and [Cuffdiff](https://cole-trapnell-lab.github.io/cufflinks/cuffdiff/)
- [Trinity2.1.1 source](https://github.com/trinityrnaseq/trinityrnaseq/blob/v2.1.1/Trinity)
- [RSEM reference preparation](https://deweylab.github.io/RSEM/rsem-prepare-reference.html), [quantification](https://deweylab.github.io/RSEM/rsem-calculate-expression.html) and [Trinity gene-map utility](https://github.com/deweylab/RSEM#generate-transcript-to-gene-map-from-trinity-output)
- [PRINSEQ parameters](https://prinseq.sourceforge.net/manual.html)
