# Figure Explanations And Validation

## Figure 1: DTE Candidate Counts

File:

`figure_1_dte_candidate_counts.svg`

What it shows:

This figure counts protein-coding transcript isoforms whose total RNA abundance
changed after UV.

DTE means:

```text
Differential Transcript Expression
```

The y-axis is candidate count. The x-axis is recovery time after UV:

```text
12 min, 30 min, 60 min
```

`up` means the isoform has higher RNA abundance in UV than in non-UV control.
`down` means the isoform has lower RNA abundance in UV than in non-UV control.

Validation:

The numbers are copied from:

`preliminary_data/dte_protein_coding_candidate_counts.csv`

and independently match:

`results/statistical_DTE_candidate_protein_coding_transcripts.csv`

Counts:

```text
12 min: 1113 down, 64 up
30 min: 574 down, 98 up
60 min: 879 down, 237 up
```

Interpretation:

UV exposure is associated with many more protein-coding isoform expression
decreases than increases, especially at 12 and 60 minutes.

## Figure 2: DTU Candidate Counts

File:

`figure_2_dtu_candidate_counts.svg`

What it shows:

This figure counts protein-coding transcript isoforms whose usage fraction
inside their gene changed after UV.

DTU means:

```text
Differential Transcript Usage
```

Usage means:

```text
isoform TPM / total TPM of all isoforms from the same gene
```

`usage_up` means the isoform takes a larger share of its gene's transcript pool
after UV.

`usage_down` means the isoform takes a smaller share of its gene's transcript
pool after UV.

Validation:

The numbers are copied from:

`preliminary_data/dtu_protein_coding_candidate_counts.csv`

and independently match:

`results/statistical_DTU_candidate_protein_coding_transcripts.csv`

Counts:

```text
12 min: 554 usage_down, 542 usage_up
30 min: 518 usage_down, 562 usage_up
60 min: 632 usage_down, 746 usage_up
```

Interpretation:

DTU is more balanced than DTE. This suggests UV may not only turn isoforms up or
down globally, but may also shift isoform preference within genes.

## Figure 3: Top DTU Usage Changes

File:

`figure_3_top_dtu_usage_changes.svg`

What it shows:

This figure shows the strongest protein-coding DTU candidates by absolute
change in usage fraction.

Each row is one transcript isoform. The horizontal axis runs from 0 to 1:

```text
0 = isoform is not used within that gene
1 = isoform accounts for almost all expression from that gene
```

The blue point is mean control usage. The UV point shows mean usage after UV.
The `delta_usage` label is:

```text
UV mean usage - control mean usage
```

Example:

```text
EGR2-201 at 30 min:
control usage = 0.00
UV usage = 0.91
delta_usage = +0.91
```

That means EGR2-201 goes from almost absent within the EGR2 transcript pool to
being the dominant EGR2 isoform after UV at 30 minutes.

Validation:

The plotted rows come from:

`preliminary_data/top_dtu_protein_coding_candidates.csv`

which is sorted from:

`results/statistical_DTU_candidate_protein_coding_transcripts.csv`

by largest absolute `delta_usage`.

## Important Scientific Caution

These figures show preliminary candidates, not final publication-level
significant isoform switches.

Reason:

With 3 control and 3 UV replicates, the exact permutation test has only 20
possible label arrangements. This makes p-value and FDR resolution limited.

The figures are good for preliminary data and hypothesis generation. For final
analysis, we should validate the strongest findings with a specialist isoform
tool such as DRIMSeq, DEXSeq, satuRn, or IsoformSwitchAnalyzeR.
