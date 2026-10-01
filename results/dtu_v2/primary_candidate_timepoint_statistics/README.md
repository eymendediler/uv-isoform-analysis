# Primary candidate time-point statistics

These figures use all 12 Salmon samples (three biological samples at each of
0, 12, 30 and 60 minutes). Statistical models retain the individual samples;
means are drawn only as descriptive visual guides.

## Files

- `primary_candidate_timepoint_heatmap.png/pdf`: rows are transcripts, columns
  are time-vs-0 contrasts, cell values are mean changes in within-gene TPM
  usage, red is increase, blue is decrease, and `*` indicates global
  Benjamini-Hochberg false discovery rate adjusted q <= 0.05.
- `<GENE>_dtu_statistics.png`: left panel shows total gene TPM; right panel
  shows transcript TPM divided by total gene TPM. Faint points are individual
  samples and connected opaque points are display means.
- `primary_candidate_replicate_values.tsv`: source data for every plotted
  replicate point.
- `primary_candidate_timepoint_statistics.tsv`: effect sizes, raw likelihood-
  ratio-test p-values and global Benjamini-Hochberg adjusted q-values.

Detailed statistical interpretation is in
`../../../docs/figure_and_statistics_guide.md`.

## Candidate figures

| Gene | Replicate-aware figure | Time-localized global-BH result |
|---|---|---|
| RUNX1 | `RUNX1_dtu_statistics.png` | RUNX1b decreases and RUNX1a increases at 60 min |
| ERN1 | `ERN1_dtu_statistics.png` | reciprocal ERN1-201/208 change at 60 min |
| FER | `FER_dtu_statistics.png` | FER-201 decreases at 12/30/60 min; shorter transcripts increase at selected times |
| EPC1 | `EPC1_dtu_statistics.png` | reciprocal EPC1-202/205 change at 30/60 min |
| SFT2D2 | `SFT2D2_dtu_statistics.png` | omnibus candidate; no individual global-BH time contrast |
| E2F3 | `E2F3_dtu_statistics.png` | reciprocal E2F3-201/202 change at 30/60 min |
| ELK4 | `ELK4_dtu_statistics.png` | omnibus candidate; no individual global-BH time contrast |
| PDK3 | `PDK3_dtu_statistics.png` | reciprocal PDK3-201/203 change at 12/30 min |
| MECP2 | `MECP2_dtu_statistics.png` | reciprocal MECP2-201/216 change at 12/30 min |
| NUFIP2 | `NUFIP2_dtu_statistics.png` | reciprocal NUFIP2-201/202 change at 60 min |
| FGF2 | `FGF2_dtu_statistics.png` | omnibus candidate; no individual global-BH time contrast |
| IER3 | `IER3_dtu_statistics.png` | omnibus candidate; no individual global-BH time contrast |

“Reciprocal” refers to within-gene transcript usage: when one transcript's
share increases, another transcript's share may decrease. It does not by itself
imply that total gene abundance is unchanged.

Biological interpretation and literature evidence are documented separately in
`../../../docs/biological_interpretation.md`. Absence of direct isoform
literature is reported as a knowledge gap rather than converted into evidence.
