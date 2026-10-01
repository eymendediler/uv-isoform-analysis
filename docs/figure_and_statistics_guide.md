# DTU figures and statistical tests: how to read and explain them

## What data are shown?

The figures start from the 12 deposited Salmon `quant.sf` files: four time
points (0, 12, 30 and 60 minutes) with three biological samples per time point.
For every sample and transcript, Salmon TPM (Transcripts Per Million) is divided
by the sum of TPM values for all annotated transcripts of that gene:

```text
within-gene usage = transcript TPM / sum of all transcript TPMs in the gene
```

This creates 12 separate usage observations for each transcript. The statistical
model is not fitted to four averages.

## How to read an individual gene figure

- Left panel, black points: total gene TPM in each of the three samples.
- Left panel, blue point/line: arithmetic mean of the three displayed TPM
  values at a time point. This line is descriptive only.
- Right panel, faint coloured points: the transcript's three separate
  within-gene usage observations at that time point.
- Right panel, coloured point/line: arithmetic mean of the three displayed
  usage observations. Again, this is a visual summary, not the model input.
- Each colour: one robust protein-coding transcript.
- `omnibus stageR gene q`: evidence that at least one within-gene transcript
  usage differs somewhere across 0, 12, 30 and 60 minutes.
- `*`: the indicated time versus 0-minute DRIMSeq contrast passed a global
  Benjamini–Hochberg (BH) false discovery rate threshold of q <= 0.05.

The mean line is present only because twelve partly overlapping points are hard
to follow visually. Removing the line would not change any p-value or q-value.

## How to read the heatmap

- One row: one GENCODE transcript from a primary candidate gene.
- Columns: 12-vs-0, 30-vs-0 and 60-vs-0 minutes.
- Number in a cell: mean TPM-usage change, `mean(time) - mean(0)`. It is an
  effect-size summary. `+0.25` means a 25 percentage-point increase in the
  transcript's share of gene TPM.
- Red: increased within-gene usage relative to 0 minutes.
- Blue: decreased within-gene usage relative to 0 minutes.
- White: little descriptive change.
- `*`: formal time-specific DRIMSeq test has global BH q <= 0.05.
- No star does not prove no change. It means the data did not localize that
  particular time contrast at the selected genome-wide error-control level.

## Why compare time differences with replicate variability?

There are two kinds of differences in the observed values:

1. **Within-time biological variation:** three independently prepared samples
   at the same time are not numerically identical.
2. **Between-time signal:** usage may systematically shift after UV.

The relevant question is not merely whether two means differ. It is whether the
observed separation between time groups is larger than would be plausible from
the amount of sample-to-sample variation seen in the experiment. Three tightly
clustered values shifting together provide stronger evidence than the same
mean shift produced by highly scattered, overlapping values.

Means are therefore allowed as descriptive summaries. What is avoided is
replacing the three observations by one mean before fitting the statistical
model, because that would discard the information needed to estimate biological
variability.

## What does Dirichlet–multinomial mean?

For a gene with several transcripts, the data for one sample are a vector of
counts, for example `[A=60, B=30, C=10]`, whose components compete for the same
gene total. A multinomial distribution models this compositional allocation.
Real biological samples vary more than a simple multinomial expects. The
Dirichlet component allows the underlying A/B/C proportions themselves to vary
between samples. The combined Dirichlet–multinomial model therefore represents:

- counts that must be allocated among transcripts of the same gene; and
- extra biological variability between replicate samples.

DRIMSeq fits regression coefficients for time while estimating a precision
parameter: lower precision corresponds to greater extra sample-to-sample
variation. This is a statistical generalized linear model, not a machine
learning classifier.

## Omnibus and time-specific tests

The primary full model is `~ time`; the null model is `~ 1`.

- `~ time` permits a different composition at 0, 12, 30 and 60 minutes.
- `~ 1` is the no-time-effect hypothesis: one common expected composition.
- The omnibus likelihood-ratio test (LRT) asks whether allowing time-specific
  compositions improves the fit enough to reject no time effect.

This ordering is appropriate for the hypothesis because the response need not
be monotonic: an isoform can rise at 12 minutes and return by 60 minutes. Testing
all four times first avoids assuming a particular temporal shape. Only after
that primary screen do the secondary 12-vs-0, 30-vs-0 and 60-vs-0 contrasts
localize the pattern.

## Where do p-values come from?

DRIMSeq calculates a likelihood-ratio statistic:

```text
2 × [log-likelihood of the more flexible model
     - log-likelihood of the restricted model]
```

Under the null hypothesis, this statistic is compared with its reference
chi-square distribution. The p-value is the probability, under the no-effect
hypothesis and model assumptions, of obtaining a likelihood-ratio statistic at
least as large as the observed one. We do not manually choose p-values. We
choose an error-control threshold before interpreting them.

## What is Benjamini–Hochberg correction?

BH means **Benjamini–Hochberg false discovery rate correction**. With `m` tests,
the raw p-values are sorted from smallest to largest. A p-value at rank `i` is
compared with its rank relative to the total number of tests; adjusted values
are made monotonic. Informally, a small raw p-value must provide stronger
evidence when it is one among tens of thousands of opportunities for a chance
result.

For the secondary localization analysis, correction was performed once across
all 75,477 transcript-by-time tests—not only across the visually attractive
candidate genes. At global BH q <= 0.05, 316 transcript-time comparisons in 158
genes remained. A q-value is an FDR-adjusted p-value; it is not the probability
that this particular result is false.

## What is the proportion-SD filter?

For every transcript, its within-gene count proportion is calculated separately
in all 12 samples and the standard deviation (SD) of these 12 proportions is
calculated without using time labels. A low SD means the transcript's observed
share hardly varies anywhere in the dataset.

The official Bioconductor rnaseqDTU workflow reported that DRIMSeq transcript
tests could exceed their nominal false discovery rate and that setting p-values
to 1 for transcripts with proportion SD below 0.10 made simulated error control
more conservative. It does not make low-SD data “better”; it prevents near-flat
transcripts from being called significant because of p-value calibration
problems. The 0.10 value is therefore a literature-based statistical safeguard,
not a universal biological law, and should be accompanied by sensitivity
analysis at alternative cutoffs.

Reference: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html

