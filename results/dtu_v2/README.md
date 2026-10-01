# Current differential transcript usage results

This directory is the authoritative result location. Files without
`scaled_tpm_prop_sd_0_1` in their names are sensitivity or superseded
intermediate outputs unless explicitly described below.

## Primary robust candidate universe

- `robust_stageR_transcripts_all_filters_scaled_tpm_prop_sd_0_1.tsv.gz`:
  every transcript confirmed by stageR after the proportion-SD safeguard in
  loose, medium and strict information-filter runs (187 transcripts).
- `robust_protein_sequence_details_scaled_tpm_prop_sd_0_1.tsv.gz`:
  annotated protein-coding subset with statistical evidence, continuous usage
  effect sizes, protein identifiers, lengths and sequence hashes.
- `robust_protein_sequence_gene_summary_scaled_tpm_prop_sd_0_1.tsv`:
  gene-level summary; the `multiple_distinct_protein_sequences` column defines
  the 12 primary protein-changing candidate genes.
- `sensitivity_summary_scaled_tpm_prop_sd_0_1.tsv`:
  retained and confirmed counts for all three information filters.

## Primary omnibus runs

`loose_unblocked_scaled_tpm/`, `medium_unblocked_scaled_tpm/` and
`strict_unblocked_scaled_tpm/` contain the model design, DRIMSeq fit, gene and
transcript omnibus results, stageR tables and filter summaries. Their corrected
stageR output is `stageR_adjusted_prop_sd_0_1.tsv`.

## Secondary time localization

- `timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1.tsv.gz`: all 75,477
  transcript-by-time tests, raw p-values, global Benjamini-Hochberg q-values,
  directions and effect sizes.
- `timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1_summary.tsv`: global test
  counts.
- `primary_candidate_timepoint_statistics/`: one figure for every primary
  protein-changing candidate, a summary heatmap and exact plotting source data.

## Interpretation

- `analysis_status.md`: accepted counts, incomplete models and boundaries.
- `../../docs/current_analysis_methods.md`: complete current methods.
- `../../docs/figure_and_statistics_guide.md`: explanation of every displayed
  quantity and statistical term.
- `../../docs/biological_interpretation.md`: gene/isoform literature review with
  evidence levels.

## Do not promote intermediate alternatives

Raw-count runs, unfiltered stageR tables, the failed blocked model and the
medium-only linear-time model are retained for sensitivity/development. They do
not replace the corrected categorical-time primary analysis.
