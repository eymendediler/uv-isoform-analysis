# GitHub Content QA Report

Date: 2026-07-22

Branch checked:

```text
preliminary-isoform-analysis
```

GitHub repository:

```text
https://github.com/eymendediler/uv-isoform-analysis
```

## What Was Checked

The GitHub branch was cloned into a temporary directory and checked independently
from the local working copy.

Checks performed:

1. Required documentation, script, result, figure, and validation files exist.
2. Raw input data are not committed.
3. Large full-result matrices are not committed.
4. No committed file exceeds GitHub's 100 MB file limit.
5. DTE candidate counts match the source result table.
6. DTU candidate counts match the source result table.
7. DTE and DTU summary files match the candidate result tables.
8. Figure 1 and Figure 2 count labels match their source CSV files.
9. Figure 3 delta-usage labels match the top DTU candidate table.
10. Gene-specific usage fractions sum to 1 when total gene TPM is greater than 0.
11. Each selected gene has a usage CSV and SVG plot.
12. Each selected gene has a replicate dot plot and ORF/CDS annotation table.
13. ORF/CDS summary contains all 10 selected genes.
14. Documentation avoids overclaiming final protein-level or publication-level
    significance.

## QA Result

```text
Errors: 0
Warnings: 0
```

## Verified Candidate Counts

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

## Interpretation

The GitHub branch contains a consistent preliminary isoform analysis package:

- scripts
- explanatory documentation
- candidate result tables
- preliminary figures
- gene-specific isoform plots
- replicate-level validation plots
- ORF/CDS annotation summaries
- reproducibility checks

The repository is suitable for preliminary discussion with the supervisor.

Important scientific limitation:

These outputs are preliminary candidate-level results. Specialist transcript
usage tools should still be used for final publication-level validation.
