# Public sequence metadata and sample mapping

This reconstruction uses the study reported in [Yang et al., DOI
10.1038/s41437-018-0090-1](https://doi.org/10.1038/s41437-018-0090-1).
The raw-read study is [SRP142554](https://www.ebi.ac.uk/ena/browser/view/SRP142554)
and BioProject [PRJNA451237](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA451237).
The repository contains metadata, **not sequence reads**.

## Verified archive snapshot

`config/runs.tsv` was obtained from the official ENA `filereport` endpoint on
2026-09-24. It contains 31 runs, 47 compressed FASTQ references and
145,909,488,473 advertised compressed bytes (about 145.9 GB / 135.9 GiB):

| Assay | Libraries/runs | Layout | FASTQ files |
| --- | ---: | --- | ---: |
| mRNA | 16 | PAIRED | 32 |
| miRNA | 15 | SINGLE | 15 |

The Methods report 16 mRNA tissue samples and 15 miRNA tissue samples; M2-EA
is excluded from miRNA profiling. The archive counts agree. ENA file size and
MD5 values are **archive-reported expectations**, not claims that we downloaded
and validated the reads.

The exact request URL, UTC retrieval time and SHA-256 of the source report
bytes are repeated on each manifest row. The checked-in snapshot's source
report SHA-256 is
`1b475aba56f0fd3545d49ca2f0dda94766a2ef8aee15343d799c914f47948aa2`.
Archive reports may change; a refreshed report is a new snapshot, and its
ordering alone can change the raw-report hash.

## Retrieve metadata only

Run from the repository root using Python 3.10 or newer; no third-party package
is required for this script:

```bash
python scripts/fetch_run_manifest.py --study SRP142554 --expect-runs 31 --output config/runs.refreshed.tsv
```

This command contacts only the ENA metadata endpoint. It **does not download
FASTQ files**, launch SRA Toolkit, execute a pipeline, or infer sample conditions.
It refuses to overwrite an existing manifest unless `--overwrite` is explicit.
All validation completes before the output is atomically published; failures
leave an existing manifest unchanged and clean up the temporary output.
The output directory must already exist. An offline import is also available:

```bash
python scripts/fetch_run_manifest.py --input-report saved-ena-report.tsv --study SRP142554 --expect-runs 31 --output config/runs.imported.tsv
```

An offline report must contain all columns requested by the script's `FIELDS`
constant. It is validated identically, but `retrieved_utc` is left blank because
the original acquisition time cannot be inferred from a local file. The
`source_url` describes the expected ENA query; it is not proof that a supplied
offline report actually came from ENA.

The implementation follows the official [ENA file-report
documentation](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)
and [ENA read-run field definitions](https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run).
It queries `result=read_run` using fixed requested fields and accepts study or
BioProject accessions, not arbitrary endpoint URLs. HTTPS redirects outside the
ENA report endpoint are rejected. Report size is capped at 16 MiB.

## Manifest schema and file-list handling

There is one row per run. Run, experiment, BioSample, library name, experiment
alias, reported layout/strategy, source/selection, instrument and read/base
counts remain separate columns. The original `fastq_ftp`, `fastq_md5` and
`fastq_bytes` fields are retained verbatim.

- Their semicolon-separated entries are positional: entry `i` in each field
  describes the same compressed file. Every entry is retained, including
  reports with more than two files; the parser never silently truncates lists.
- Empty lists are permitted only when **all three** fields are empty, producing
  `fastq_file_count=0`. This means no ENA-generated FASTQ is listed, not an empty
  biological sample. A listed file must have a 32-hex MD5 and positive integer
  byte size. Missing elements and mismatched list lengths fail validation.
- `fastq_https` has the same order and cardinality as `fastq_ftp`. Only locations
  on exactly `ftp.sra.ebi.ac.uk` under `/vol1/fastq/` are accepted. Scheme-less
  ENA locations and FTP locations become HTTPS; credentials, ports, other
  hosts, path traversal, URL parameters and fragments are rejected.
- `fastq_file_count` and `fastq_total_bytes` are derived metadata. File count
  does **not** determine `library_layout`. In general, paired runs can include
  an unpaired-read file; do not assign mates by list position alone.
- Missing required columns, malformed records, duplicate run IDs, invalid
  accessions, unexpected study membership and empty reports fail closed.
  `--expect-runs 31` additionally guards the expected study run count.

## Sample-to-run map

These mappings use agreement between ENA `experiment_alias` and `library_name`,
the reported assay/layout, and the **explicit sample labels in the paper's
Methods**. The `mRNA_` prefix is removed only under those checks. The fetcher
does not use a heuristic based only on a library name. Unrecognized or
conflicting metadata leaves `sample_id` and `assay` blank with
`mapping_evidence=unresolved`; do not use those rows without review.

| Paper sample | mRNA run | miRNA run |
| --- | --- | --- |
| F1-O | SRR7062132 | SRR7062137 |
| F2-O | SRR7062129 | SRR7062144 |
| F3-O | SRR7062130 | SRR7062143 |
| F4-O | SRR7062127 | SRR7062153 |
| F5-O | SRR7062128 | SRR7062152 |
| F6-O | SRR7062125 | SRR7062151 |
| F2-E | SRR7062126 | SRR7062150 |
| F3-E | SRR7062145 | SRR7062149 |
| M1-OA | SRR7062148 | SRR7062140 |
| M1-OB | SRR7062147 | SRR7062139 |
| M2-OA | SRR7062146 | SRR7062142 |
| M2-OB | SRR7062155 | SRR7062141 |
| M3-OA | SRR7062154 | SRR7062136 |
| M3-OB | SRR7062133 | SRR7062135 |
| M2-EA | SRR7062134 | Not included in the paper |
| M2-EB | SRR7062131 | SRR7062138 |

**Archive caveat:** all 31 runs share BioSample `SAMN08965226` / secondary
sample accession `SRS3217376`. The same multi-sample `sample_alias` occurs across
the deposit. Neither this BioSample accession nor that alias identifies an
individual animal or a tissue sample. Do not group biological replicates by
them. `config/samples.tsv` holds the separately curated experimental design.

Condition, tissue, species and individual identity cannot be established from
library layout or library names alone. The paper explains the design:
Finnsheep ovarian samples come from six animals; mouflon OA/OB are samples from
only three animals, not six independent animals. Mouflon EA/EB are from the
two uterine horns of **one** animal, M2. F2/F3 also contribute samples across
tissues. Finnsheep ovarian sampling was at the follicular growth phase and
endometrial sampling in early pregnancy; mouflon cycle stage was unknown
except that M2 had a corpus luteum. The design does not establish matched
reproductive stages between species.

## Deliberate read staging

Downloading about 146 GB of compressed data requires a separate, deliberate
action and additional space for intermediate results. No download mode is
included here. For any independently downloaded archive file:

1. Use the exact trusted HTTPS URL from the manifest and download to a temporary
   filename, never directly to a workflow input filename.
2. Check the **compressed file's** byte size and MD5 against the aligned manifest
   entries; an incomplete or mismatched file must not be used.
3. Verify the read roles from the run metadata/filenames and inspect FASTQ
   records; only then atomically rename or link validated files to the paths
   in `config/samples.tsv`. Record the run-to-local-path association.

The supplied sample sheet expects `data/mrna/SAMPLE_R1.fastq.gz`,
`data/mrna/SAMPLE_R2.fastq.gz` and `data/mirna/SAMPLE.fastq.gz` (with `.` for the
absent M2-EA miRNA library). It contains no downloaded files. Never substitute
one assay for another, invent an M2-EA miRNA run, or silently discard extra
archive files. The original paper also used different historical reference
assemblies for mRNA and miRNA; see the workflow's explicit reference settings.
