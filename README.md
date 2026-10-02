# UV time-course transcript isoform analysis

This repository contains a reanalysis of the HeLa-S3 RNA-seq time course from
Kaya and Adebali (2025), *UV-induced reorganization of 3D genome mediates DNA
damage response* (GEO: `GSE268349`). The biological question is whether UV
exposure is followed by time-dependent changes in within-gene transcript usage,
particularly changes between transcripts capable of encoding different protein
sequences.

## Analysis status

**Use the `dtu_v2` analysis for current biological interpretation.** Earlier
permutation and SUPPA2 outputs are retained for provenance and method
development, but are superseded and must not be presented as the primary result.

## Start here: current figures

The figure below and the linked gene panels are the current, presentation-ready
results. Files located directly under the older `preliminary_data/` and
top-level `results/` locations are retained only as legacy outputs.

![Primary time-localized DTU candidates](results/dtu_v2/primary_candidate_timepoint_statistics/primary_candidate_timepoint_heatmap.png)

- [RUNX1 replicate-aware DTU figure](results/dtu_v2/primary_candidate_timepoint_statistics/RUNX1_dtu_statistics.png)
- [E2F3 replicate-aware DTU figure](results/dtu_v2/primary_candidate_timepoint_statistics/E2F3_dtu_statistics.png)
- [FER replicate-aware DTU figure](results/dtu_v2/primary_candidate_timepoint_statistics/FER_dtu_statistics.png)
- [EPC1 replicate-aware DTU figure](results/dtu_v2/primary_candidate_timepoint_statistics/EPC1_dtu_statistics.png)
- [All 12 current candidate figures and source tables](results/dtu_v2/primary_candidate_timepoint_statistics/)

- [Current methods](docs/current_analysis_methods.md)
- [Exact reproducibility commands](docs/reproducibility_commands.md)
- [Current analysis status and accepted result counts](results/dtu_v2/analysis_status.md)
- [Current figures and source data](results/dtu_v2/primary_candidate_timepoint_statistics/README.md)
- [Figure and statistical interpretation guide](docs/figure_and_statistics_guide.md)
- [Candidate-by-candidate biological interpretation](docs/biological_interpretation.md)
- [Questions to confirm with the original laboratory](docs/questions_for_supervisor.md)
- [Superseded/legacy analysis notice](LEGACY_ANALYSES.md)

## Experimental data

The analysis retains all 12 samples:

| Recovery after UV | Samples |
|---|---:|
| 0 min / non-UV control | 3 |
| 12 min | 3 |
| 30 min | 3 |
| 60 min | 3 |

Inputs are the deposited Salmon transcript quantifications and GENCODE v35
annotation. The current model treats samples as independent because the public
metadata do not establish that replicate labels form matched experimental
blocks across time. A blocked sensitivity analysis should be attempted only if
the original laboratory confirms that relationship.

## Why an omnibus time model?

The temporal form of an isoform response is not known in advance. Biologically
plausible patterns include:

- a transient response at 12 minutes;
- a response beginning at 30 minutes and persisting at 60 minutes;
- a monotonic increase or decrease; and
- an early increase followed by recovery.

The primary categorical-time omnibus test therefore asks whether any of the
four time points differ without imposing a linear or monotonic trajectory.
Secondary `12 vs 0`, `30 vs 0`, and `60 vs 0` contrasts localize the supported
time points.

## Current statistical workflow

```text
12 Salmon quantifications
  -> tximport length-scaled count-like values (countsFromAbundance="scaledTPM")
  -> DRIMSeq information filtering
  -> DRIMSeq Dirichlet-multinomial model: full ~ time, null ~ 1
  -> omnibus likelihood-ratio test
  -> condition-blind transcript proportion-SD filter
  -> stageR gene screening and transcript confirmation (target OFDR = 0.05)
  -> robustness intersection across loose, medium and strict filters
  -> GENCODE protein-sequence comparison
  -> secondary time-vs-0 DRIMSeq contrasts
  -> global Benjamini-Hochberg correction across all transcript x time tests
```

The current corrected primary analysis contains 187 stageR-confirmed
transcripts shared by all three information-filter settings; 111 are
protein-coding. Twelve genes contain at least two robust coding transcripts
with distinct GENCODE amino-acid sequences. These are RNA-level candidates,
not proof of translated protein switching.

## Primary candidate figures

Replicate-aware plots and source tables are available for:

`RUNX1`, `ERN1`, `FER`, `EPC1`, `SFT2D2`, `E2F3`, `ELK4`, `PDK3`, `MECP2`,
`NUFIP2`, `FGF2`, and `IER3`.

In all plots, individual points are the biological samples. Connected means are
descriptive visual guides only; the statistical model is fitted to individual
sample counts and does not replace replicates with means.

## Reproducing the current analysis

The versioned environment is described in `environment-dtu.yml`. Input paths
and checksums are documented under `metadata/`.

```bash
conda env create --prefix .envs/dtu -f environment-dtu.yml
conda run --prefix .envs/dtu Rscript scripts/run_drimseq_omnibus.R \
  --filter loose --model unblocked --input-scale scaled_tpm
conda run --prefix .envs/dtu Rscript scripts/run_drimseq_omnibus.R \
  --filter medium --model unblocked --input-scale scaled_tpm
conda run --prefix .envs/dtu Rscript scripts/run_drimseq_omnibus.R \
  --filter strict --model unblocked --input-scale scaled_tpm
conda run --prefix .envs/dtu Rscript scripts/apply_drimseq_proportion_sd_filter.R
python scripts/summarize_dtu_sensitivity.py
python scripts/validate_protein_sequences.py
conda run --prefix .envs/dtu Rscript scripts/run_drimseq_timepoint_contrasts.R
python scripts/plot_primary_dtu_statistics.py
```

Alternatively, run the complete accepted sequence with:

```bash
bash scripts/run_current_analysis.sh .envs/dtu
```

See [the exact command record](docs/reproducibility_commands.md) and
[scripts/README.md](scripts/README.md) for the role and status of every current
and superseded script.

## Reference study

Kaya VO, Adebali O. UV-induced reorganization of 3D genome mediates DNA damage
response. *Nature Communications* 16, 1376 (2025).
https://www.nature.com/articles/s41467-024-55724-7
