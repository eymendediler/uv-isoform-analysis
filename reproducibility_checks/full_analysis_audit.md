# Full analysis audit

Date: 2026-07-28

This check covered the sample map, raw Salmon files, DTE/DTU calculations,
SUPPA2 inputs and outputs, summary tables, and replicate-level plots.

## Samples

All 12 expected Salmon files are present:

| Condition | Replicates |
|---|---|
| control | 1, 2, 3 |
| UV 12 min | 1, 2, 3 |
| UV 30 min | 1, 2, 3 |
| UV 60 min | 1, 2, 3 |

Each file contains the same 229,580 unique transcript IDs in the same order.
There are no duplicated transcript IDs, missing values, negative TPM values or
negative estimated read counts. TPM sums are approximately one million in
every sample.

Transcript `Length` is identical in all files. `EffectiveLength` differs
slightly between samples, as expected from Salmon's sample-specific fragment
length correction. The analysis uses the first file's effective length only as
an output annotation; DTE, DTU and SUPPA2 calculations use TPM values and are
not affected by this choice.

All 229,580 quantified transcript IDs map to one GENCODE v35 transcript record.

## DTE and DTU

The candidate sets were recalculated directly from the raw Salmon TPM files
without reading the saved candidate flags.

- DTE candidate transcript sets match the saved tables exactly at 12, 30 and
  60 min.
- DTU candidate transcript sets match the saved tables exactly at all three
  time points.
- The protein-coding subsets also match exactly.
- DTE uses three control and three UV TPM columns in every comparison.
- DTU uses three control and three UV usage columns in every comparison.
- Per-sample isoform usage sums to 1 for every expressed annotated gene, within
  floating-point precision.
- Exact permutation candidate p-values are 0.10, which is the minimum
  two-sided value for three-versus-three data with 20 label assignments.

Saved candidate totals:

| Analysis | All candidates | Protein-coding candidates |
|---|---:|---:|
| DTE | 7,115 | 2,965 |
| DTU | 5,651 | 3,554 |

All DTE/DTU summary CSV files reproduce the counts in their source candidate
tables.

## SUPPA2

The four TPM condition matrices and four transcript-PSI matrices each contain
exactly three replicate columns with matching sample names. Transcript IDs are
unique.

All three transcript-DTU comparisons and all 21 local-event comparisons
(seven event types at three time points) were checked against their raw
SUPPA2 `.dpsi` files:

- candidate extraction matches the saved tables exactly;
- BH-adjusted p-values match `statsmodels` independently;
- candidate summary counts match the candidate tables;
- method-overlap and Jaccard values match an independent recalculation;
- all 4,778 overlapping DTU rows have the same direction in both methods;
- the set of 2,124 recurrent transcripts is reproduced exactly.

The entire SUPPA2 workflow and the cross-method comparison were then run three
times from the scripts. All six final SUPPA2 CSV files had identical SHA-256
checksums in the starting output and each of the three runs. This includes the
transcript-DTU candidates, local-event candidates, summaries, overlap table and
recurrent-transcript table. The checksums are recorded in
`reproducibility_checks/reproducibility_report.md`.

SUPPA2 transcript rows passing both `padj_BH <= 0.05` and
`absolute delta PSI >= 0.10`:

| Time | Rows |
|---|---:|
| 12 min | 23 |
| 30 min | 23 |
| 60 min | 46 |

Across the local-event comparisons, 32 event/time-point rows pass the same FDR
and effect-size criteria.

SUPPA2 assigns an empirical p-value of exactly zero to some of the strongest
events. Their BH-adjusted values are therefore also zero. This is an output of
SUPPA2's finite empirical null procedure, not evidence that the biological
false-positive probability is literally zero.

## Plots and preliminary tables

- DTE and DTU preliminary count tables match their source result tables.
- All ten selected genes have a usage CSV, usage SVG, replicate SVG and
  ORF/CDS table.
- Usage fractions sum to 1 whenever the gene has nonzero expression.
- Gene-specific stacked bars include all transcript expression. Transcripts
  below the individual display threshold are combined as `Other expressed
  isoforms`, so every expressed condition reaches a total usage fraction of 1.
- Every replicate plot contains three points per condition for each displayed
  isoform (12 points per isoform across four conditions).
- The publication figures were regenerated directly from the saved result
  tables and sample-level TPM values. Their source-data CSV files match those
  inputs, and every gene/transcript/condition group in Figure 5 contains
  replicates 1, 2 and 3 exactly once.
- Publication figures are available as editable SVG, single-page PDF and
  300-dpi PNG files. The exported files were checked for successful rendering.
- All Python scripts compile successfully.

## Audit conclusion

The sample assignment, use of three replicates per condition, numerical
calculations and saved result tables are internally consistent. No result-table
error was found.

One documentation error was found and corrected: earlier text stated that no
SUPPA2 candidates passed BH correction. The correct counts are 23, 23 and 46
transcript/time-point rows after applying both the FDR and delta-PSI cutoffs.

The statistical limitations remain unchanged: each comparison has only three
replicates per condition, transcript abundance is estimated from short reads,
the custom permutation test has low p-value resolution, and SUPPA2's empirical
zero p-values should be interpreted cautiously.
