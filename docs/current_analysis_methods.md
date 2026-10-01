# Current analysis methods

## Objective

The primary hypothesis is that UV exposure is followed by time-dependent
changes in transcript composition within genes in HeLa-S3 cells. The analysis
prioritizes genes for which multiple supported transcripts encode distinct
amino-acid sequences, while recognizing that short-read RNA sequencing cannot
establish protein production.

## Samples and design

Twelve Salmon quantifications were retained: three samples each at 0, 12, 30
and 60 minutes after UV exposure. Time is represented as a categorical factor.
The public metadata do not establish that equal replicate numbers across time
are matched experimental blocks, so the accepted primary model is unblocked.

The categorical design intentionally permits transient, delayed, monotonic and
recovery patterns. No single temporal trajectory is imposed.

## Quantification import

Salmon `quant.sf` files are imported with tximport at transcript resolution
using `countsFromAbundance="scaledTPM"`. These length-scaled count-like values
preserve a count scale suitable for modeling while reducing transcript-length
bias in within-gene composition. Salmon TPM is used separately for descriptive
within-gene usage effect sizes.

## Information filtering

Loose, medium and strict DRIMSeq filters were run independently. They differ in
minimum transcript count and proportion requirements. A result is called robust
only if it is confirmed under all three filter settings. Filter choices and
retained feature counts are reported rather than silently selecting one
threshold after viewing candidates.

## Dirichlet-multinomial model

For one sample, the counts of all transcripts of a gene form a composition:
increasing the share allocated to one transcript changes the shares available
to the others. A multinomial distribution models this allocation conditional
on the gene total. A simple multinomial assumes one fixed composition for all
replicates. The Dirichlet layer instead treats each sample's underlying
composition as a draw around a group-level composition, allowing biological
overdispersion. DRIMSeq estimates both time coefficients and precision; lower
precision represents greater extra sample-to-sample variability.

This model is fitted to individual sample values. Arithmetic means shown in
figures are calculated only after the replicate-level data are retained and
are descriptive visual guides; no statistical test replaces three biological
samples with one mean.

## Primary omnibus test

The full model and null model are:

```text
full: ~ time
null: ~ 1
```

The full model allows a separate expected transcript composition at each time.
The null model represents one common expected composition. DRIMSeq performs a
likelihood-ratio test (LRT):

```text
LRT = 2 x [log-likelihood(full model) - log-likelihood(null model)]
```

The LRT is compared with its reference chi-square distribution. The p-value is
the probability, under the no-time-effect hypothesis and model assumptions, of
obtaining an LRT statistic at least as large as the observed statistic.

The omnibus test is primary because the biological response may be transient,
delayed, monotonic or recovering. It asks whether any time differs without
requiring a prespecified response shape.

## Condition-blind proportion-SD safeguard

For each retained transcript, its count proportion within its gene is
calculated in every sample. The standard deviation (SD) of these 12 proportions
is calculated without using time labels. Following the Bioconductor rnaseqDTU
workflow, transcript tests with proportion SD below 0.10 are assigned p=1
before stageR adjustment. This safeguard improved transcript-level FDR/OFDR
calibration in the workflow's simulations; it is not claimed as a universal
biological cutoff.

## Hierarchical error control

stageR first screens genes for any evidence of differential transcript usage
and then confirms transcripts within screened genes. The target overall false
discovery rate (OFDR) is 0.05. The accepted primary results are transcripts
confirmed by stageR after the proportion-SD safeguard in all three information
filter settings.

## Protein-sequence prioritization

Robust transcripts are mapped to GENCODE v35 protein identifiers. Amino-acid
sequences are compared by SHA-256 digest, so different protein identifiers that
encode identical sequences are not counted as distinct protein products. A
primary protein-changing candidate gene must have at least two robust coding
transcripts with different amino-acid sequences.

## Secondary time-point localization

The saved strict-filter DRIMSeq fit is tested with coefficients `time12`,
`time30` and `time60`, each relative to the 0-minute reference. The same
condition-blind proportion-SD safeguard is applied. Benjamini-Hochberg false
discovery rate correction is performed once across all 75,477 transcript-by-
time tests. These contrasts localize the time pattern but do not replace the
primary omnibus/stageR discovery analysis.

## Effect sizes and figures

For sample `s`, transcript `t`, and gene `g`:

```text
TPM usage(t,s) = TPM(t,s) / sum[TPM(all annotated transcripts of g,s)]
```

Heatmap cell values are `mean usage(time) - mean usage(0)` and are descriptive
effect sizes. Individual gene figures show all three biological samples at each
time. Connected means are visual summaries only. Stars identify secondary
contrasts with global Benjamini-Hochberg q <= 0.05.

## Interpretation boundary

The analysis establishes RNA-level time-associated differential transcript
usage under the stated model. It does not establish translation, protein
stability, protein localization, causal involvement in the DNA damage response,
or a complete full-length isoform structure. Those claims require orthogonal
validation, ideally transcript-specific junction assays and targeted long-read
or protein-level measurements.

## Method references

- DRIMSeq: https://pmc.ncbi.nlm.nih.gov/articles/PMC5200948/
- rnaseqDTU workflow: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html
- stageR: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1277-0
- tximport: https://f1000research.com/articles/4-1521
- RNA-seq design and blocking: https://pmc.ncbi.nlm.nih.gov/articles/PMC2881125/
