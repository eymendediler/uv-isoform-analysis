# Publication figures

## Figure 4: analysis overview

`figure_4_analysis_overview` summarizes the main candidate counts:

- protein-coding DTE candidates;
- protein-coding DTU candidates;
- overlap between the permutation DTU analysis and SUPPA2;
- SUPPA2 local splicing event classes.

The percentages in panel C use the number of common transcript IDs divided by
the candidate total from each method. The number printed near the bottom of
each time point is the absolute overlap.

## Figure 5: selected candidate genes

`figure_5_candidate_isoform_usage` shows replicate-level isoform usage for
ZNF831, RUNX1, SKIL and HDAC11. Open circles are individual biological
replicates and filled points connected by lines are condition means.

Usage is calculated separately in every sample:

```text
transcript TPM / total TPM of all transcripts from the same gene
```

The figure is intended to show both the effect size and the variability between
the three replicates. It does not replace a formal significance test.

## How to read the figures

Figure 4 gives the broad view. Panels A and B show how many protein-coding
transcripts passed the preliminary DTE and DTU filters. Panel C asks whether
the permutation and SUPPA2 approaches find the same DTU candidates. Panel D
shows which local splicing-event classes contribute most candidates.

Figure 5 is a closer look at four genes. Each open circle is one biological
replicate, so three circles should be visible for every transcript and
condition. The filled point is their mean. A rising line means that the
transcript forms a larger fraction of its gene's expression after UV; a
falling line means that its fraction becomes smaller.

Examples:

- In ZNF831, usage shifts progressively from ZNF831-201 toward ZNF831-203.
  The change is already visible at 12 min and is strongest at 60 min.
- RUNX1-202 decreases after UV, while RUNX1-203 increases, particularly at
  60 min. This is consistent with a change in relative isoform usage.
- SKIL-201 decreases across the time course, whereas SKIL-205 increases and
  then remains near the same level at 30 and 60 min.
- HDAC11 shows a transient pattern: HDAC11-213 is prominent at 30 min, while
  HDAC11-201 is prominent in control, 12 min and 60 min samples.

These plots support candidate-level observations, not a claim of a confirmed
biological mechanism. The small sample size, short-read transcript estimates
and statistical thresholds should be considered when interpreting them.

## Formats and source data

Each figure is provided as:

- SVG for editable vector graphics;
- PDF for manuscript or presentation use;
- 300-dpi PNG for quick viewing.

The CSV files in this directory contain the exact values used to draw each
panel. Figures can be regenerated with:

```bash
python scripts/create_publication_figures.py
```
