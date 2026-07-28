# Isoform Analysis Reproducibility Check

Date: 2026-07-22

The isoform analysis pipeline was run three times:

1. Existing/reference run
2. Repeat run 2
3. Repeat run 3

The following scripts were rerun for repeat runs 2 and 3:

```text
scripts/statistical_isoform_analysis.py
scripts/create_preliminary_isoform_report.py
scripts/create_gene_specific_isoform_plots.py
scripts/create_isoform_validation_outputs.py
```

## Main Candidate Counts

Protein-coding DTE candidates:

```text
12 min: 1113 down, 64 up
30 min: 574 down, 98 up
60 min: 879 down, 237 up
```

Protein-coding DTU candidates:

```text
12 min: 554 usage_down, 542 usage_up
30 min: 518 usage_down, 562 usage_up
60 min: 632 usage_down, 746 usage_up
```

## Checksum Validation

These SHA-256 checksums were identical across the reference run, repeat run 2,
and repeat run 3.

```text
e95f7f3b83ba887b46289bc501e0204e47de3497457ebb67692552aa962b88df  results/statistical_DTE_candidate_protein_coding_transcripts.csv
3472c356e2ec7831a4cfb56937830ccb30b4e189d1500e37eee0c0389a724e74  results/statistical_DTU_candidate_protein_coding_transcripts.csv
00f5c85e4ee1d95a1e65703dc66e5eb1e38388a9be6cb1b21aabd111d363b8ca  preliminary_data/dte_protein_coding_candidate_counts.csv
dc4053da488ea966e780c734eaaa4ea5852cb40b694b54a1126b8ed4568c87c2  preliminary_data/dtu_protein_coding_candidate_counts.csv
6343188231910e321dd2bc5df191cc689b3439d776ae7d322ade7f10686b88c0  preliminary_data/isoform_validation/selected_genes_orf_length_summary.csv
```

## Interpretation

The pipeline is reproducible for the current input files and code. Re-running
the analysis produces identical candidate tables, summary counts, and ORF/CDS
summary outputs.

Important limitation:

These are preliminary candidate-level results. The exact permutation test is
limited by the 3 control and 3 UV replicates per comparison. Final validation
should use specialist transcript usage tools such as DRIMSeq, DEXSeq, satuRn, or
IsoformSwitchAnalyzeR.

## SUPPA2 Repeat Check

Date: 2026-07-28

After the SUPPA2 analysis was added, the complete SUPPA2 workflow was run three
times. Every run rebuilt the transcript TPM and PSI matrices, transcript-level
DTU comparisons, seven local splicing event classes at three time points, and
the cross-method comparison tables.

The following SHA-256 checksums were identical in the starting output and all
three repeat runs:

```text
d4f049cfb03d5955191344569b259b66b33f0690a236b99244559490312a4cbf  results/suppa2/cross_method_supported_isoforms.csv
81a25289d0bfe207a66fe1bf5596093b16cb6f012e4a12fb9b0922a8c3f34588  results/suppa2/method_concordance_summary.csv
1730da12f754e18a7cb982235c72526662e0b7ed07ca07d8839408f81c4e5543  results/suppa2/recurrent_isoform_dtu_candidates.csv
270b69e92d222ec8eb643c97665f68f236a23f865715fbfbf5f472fc6af98112  results/suppa2/suppa2_candidate_summary.csv
96a50644e15eb4614bb8e650fed597c17e25022b36439008d62647df1208af51  results/suppa2/suppa2_isoform_dtu_candidates.csv
708bab2af375433fefc11bb7ebcbca9447741dcc409fe90a6d5805f8f08d66c1  results/suppa2/suppa2_local_splicing_candidates.csv
```

This confirms that the SUPPA2 and cross-method result generation is
deterministic for the current inputs, code and tool version.

The run-by-run manifest is available at:

```text
reproducibility_checks/suppa2_repeat_checksums.csv
```
