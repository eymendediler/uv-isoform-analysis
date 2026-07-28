# SUPPA2 analysis

I used SUPPA2 v2.4 as an independent check of the permutation-based DTU
analysis. The input consists of Salmon TPM estimates and the GENCODE v35
annotation.

Three comparisons were run:

```text
control vs 12 min
control vs 30 min
control vs 60 min
```

SUPPA2 was used at two levels:

- transcript PSI, for differential transcript usage;
- local splicing PSI for skipped exon (SE), alternative splice sites (A3/A5),
  mutually exclusive exon (MX), retained intron (RI), and alternative
  first/last exon events (AF/AL).

## Filters

The working candidate cutoff was:

```text
absolute delta PSI >= 0.10
SUPPA2 empirical p-value <= 0.05
average TPM >= 1
```

The tables also contain BH-adjusted p-values. After requiring both
`padj_BH <= 0.05` and `absolute delta PSI >= 0.10`, there are 23, 23 and 46
transcript/time-point rows at 12, 30 and 60 min. SUPPA2 reports an empirical
p-value of exactly zero for the strongest events, so this FDR subset should
still be treated cautiously rather than as a final validation set.

## Results

Transcript-DTU candidate counts:

| Time | Candidates |
|---|---:|
| 12 min | 2,739 |
| 30 min | 2,734 |
| 60 min | 3,387 |

Comparison with the exact-permutation DTU results:

| Time | Common transcripts | Jaccard |
|---|---:|---:|
| 12 min | 1,429 | 0.476 |
| 30 min | 1,426 | 0.473 |
| 60 min | 1,923 | 0.517 |

There were 4,778 overlapping transcript/time-point rows in total. Delta PSI
and delta usage had the same sign in every overlapping row. This agreement is
useful for ranking, although both methods start from the same Salmon estimates
and are therefore not fully independent measurements.

SUPPA2 also returned 2,124 transcripts at two or more time points. The local
analysis produced 8,190 candidate event/time-point rows. AF and SE were the
largest event classes.

## Running the analysis

```bash
python -m pip install -r requirements.txt
git clone --depth 1 https://github.com/comprna/SUPPA.git .tools/SUPPA
python scripts/run_suppa2_analysis.py
python scripts/compare_isoform_methods.py
```

The three samples in each condition are biological replicates. Repeating the
workflow is a computational reproducibility check and does not increase the
number of biological replicates.

The runner makes one performance-only change to SUPPA2: a linear lookup in a
sorted list is replaced with the equivalent binary search. It checks for the
expected v2.4 source line before making the replacement.

Large PSI, TPM and event matrices are written to `results/suppa2_work/` and are
excluded from Git.

## Output tables

```text
results/suppa2/suppa2_isoform_dtu_candidates.csv
results/suppa2/suppa2_local_splicing_candidates.csv
results/suppa2/suppa2_candidate_summary.csv
results/suppa2/method_concordance_summary.csv
results/suppa2/cross_method_supported_isoforms.csv
results/suppa2/recurrent_isoform_dtu_candidates.csv
```

For follow-up, I would start with protein-coding transcripts that are supported
by both methods, recur at more than one time point, have a large usage change,
and show consistent replicate-level expression.
