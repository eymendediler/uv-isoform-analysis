# UV-associated transcript isoform analysis

This repository contains a preliminary isoform-level reanalysis of the RNA-seq
data from Kaya and Adebali (2025), *UV-induced reorganization of 3D genome
mediates DNA damage response*.

The published study reports gene-level RNA-seq results. Here I used the Salmon
quantifications deposited under GEO accession `GSE268349` to examine transcript
expression and isoform usage after UV treatment.

## Data

- control, 12, 30 and 60 min samples (three replicates per condition)
- Salmon `quant.sf.gz` files from `GSE268349`
- GENCODE v35 transcript annotation

Raw data and the GENCODE GTF are not included in the repository.

## Analysis

I looked at three related questions:

- Differential transcript expression (DTE): whether the abundance of an
  individual transcript changes after UV.
- Differential transcript usage (DTU): whether a transcript's share of the
  total expression of its gene changes.
- SUPPA2 analysis: an independent DTU analysis together with local splicing
  events (SE, A3/A5, MX, RI and AF/AL).

The initial DTE/DTU analysis uses an exact label-permutation test for each
control-versus-UV comparison. SUPPA2 uses replicate-level PSI values and its
empirical significance procedure.

## Main results

Protein-coding candidates from the permutation analysis:

| Time | DTE down | DTE up | DTU down | DTU up |
|---|---:|---:|---:|---:|
| 12 min | 1,113 | 64 | 554 | 542 |
| 30 min | 574 | 98 | 518 | 562 |
| 60 min | 879 | 237 | 632 | 746 |

SUPPA2 identified 2,739, 2,734 and 3,387 transcript-DTU candidates at 12, 30
and 60 min, respectively. The overlap between SUPPA2 and the permutation DTU
analysis was 1,429, 1,426 and 1,923 transcripts. The direction of the usage
change agreed in all overlapping cases.

There were 2,124 SUPPA2 candidates present at two or more time points. The
local-event analysis returned 8,190 candidate event/time-point combinations,
with alternative first exon and skipped exon events as the largest groups.

These numbers are useful for candidate selection, but they should not be read
as a final list of significant isoform switches. There are only three
replicates per condition. Among the SUPPA2 candidates, 23, 23 and 46
transcript/time-point rows passed BH correction at 12, 30 and 60 min,
respectively. These require cautious interpretation because SUPPA2's empirical
procedure returned a raw p-value of zero for the strongest events.

## Repository structure

- `scripts/`: analysis and plotting scripts
- `results/`: DTE, DTU and SUPPA2 result tables
- `preliminary_data/`: summary tables and candidate plots
- `docs/`: notes on the analysis and interpretation
- `reproducibility_checks/`: rerun and content checks

The main SUPPA2 notes are in
[`docs/suppa2_isoform_analysis.md`](docs/suppa2_isoform_analysis.md).
The complete sample and result audit is in
[`reproducibility_checks/full_analysis_audit.md`](reproducibility_checks/full_analysis_audit.md).
Run-by-run SUPPA2 checksums are in
[`reproducibility_checks/suppa2_repeat_checksums.csv`](reproducibility_checks/suppa2_repeat_checksums.csv).
Publication-style overview and candidate-gene figures are in
[`preliminary_data/publication_figures/`](preliminary_data/publication_figures/).

## Reproduction

The local input files should be placed in `geo_rnaseq_quant/`, with
`gencode.v35.annotation.gtf.gz` in the repository root.

```bash
python -m pip install -r requirements.txt
python scripts/analyze_transcript_isoforms.py
python scripts/statistical_isoform_analysis.py
python scripts/create_preliminary_isoform_report.py
python scripts/create_gene_specific_isoform_plots.py
python scripts/create_isoform_validation_outputs.py
python scripts/create_publication_figures.py
```

For the SUPPA2 analysis:

```bash
git clone --depth 1 https://github.com/comprna/SUPPA.git .tools/SUPPA
python scripts/run_suppa2_analysis.py
python scripts/compare_isoform_methods.py
```

The dataset contains three biological replicates per condition. Separately,
both the preliminary DTE/DTU pipeline and the later SUPPA2 workflow were
executed in three independent computational reruns. These reruns test
reproducibility; they do not add biological replicates. Candidate tables,
summary counts, cross-method tables and ORF/CDS summaries were identical by
SHA-256 checksum.
