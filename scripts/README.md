# Analysis scripts and status

## Current primary workflow

- `audit_salmon_inputs.py`: validates sample files, transcript IDs and library
  summaries.
- `export_gencode_annotation.py`: exports the transcript-to-gene annotation.
- `run_drimseq_omnibus.R`: imports Salmon values through tximport, applies the
  named information filter, fits the categorical Dirichlet-multinomial model
  and performs the omnibus likelihood-ratio test and stageR adjustment.
- `apply_drimseq_proportion_sd_filter.R`: applies the condition-blind
  proportion standard-deviation filter recommended by the Bioconductor
  rnaseqDTU workflow, then reruns stageR.
- `summarize_dtu_sensitivity.py`: intersects results across loose, medium and
  strict filters.
- `validate_protein_sequences.py`: maps robust coding transcripts to GENCODE
  proteins and verifies distinct amino-acid sequences.
- `run_gene_expression_timecourse.R`: complementary gene-level DESeq2
  categorical-time likelihood-ratio test.
- `summarize_biological_evidence.py`: integrates statistical, effect-size and
  protein-sequence evidence.
- `run_drimseq_timepoint_contrasts.R`: secondary 12-vs-0, 30-vs-0 and 60-vs-0
  transcript contrasts with one global Benjamini-Hochberg correction family.
- `plot_primary_dtu_statistics.py`: replicate-aware candidate figures, heatmap
  and plotting source tables.

## Sensitivity/development scripts

- `dtu_filter_sensitivity.R`: independent checks of DRIMSeq filter retention.
- `plot_candidate_timecourses.py`: superseded by
  `plot_primary_dtu_statistics.py` for reporting.

## Legacy preliminary workflow

The following scripts are retained to reproduce the initial exploratory
permutation/SUPPA2 outputs but are not the current primary analysis:

- `analyze_transcript_isoforms.py`
- `statistical_isoform_analysis.py`
- `run_suppa2_analysis.py`
- `compare_isoform_methods.py`
- `create_preliminary_isoform_report.py`
- `create_gene_specific_isoform_plots.py`
- `create_isoform_validation_outputs.py`
- `create_publication_figures.py`

See `../LEGACY_ANALYSES.md` before interpreting legacy outputs.
