# DTE and DTU analysis notes

The script `scripts/statistical_isoform_analysis.py` performs separate tests for
transcript expression and transcript usage. I kept these analyses separate
because an isoform can change in absolute abundance without changing its share
of gene expression, or vice versa.

## DTE

Differential transcript expression (DTE) compares the abundance of each
transcript between control and one UV time point. TPM values are log2
transformed after adding a pseudocount of 0.1.

A transcript is retained as a DTE candidate when:

```text
maximum condition mean TPM >= 2
absolute log2 fold change >= 1
permutation p-value <= 0.10
```

## DTU

For each sample, transcript usage is calculated as:

```text
transcript TPM / sum of TPMs from the same gene
```

The reported effect size is:

```text
delta usage = mean UV usage - mean control usage
```

A transcript is retained as a DTU candidate when:

```text
maximum gene mean TPM >= 2
the gene has more than one annotated transcript
absolute delta usage >= 0.10
permutation p-value <= 0.10
```

Positive delta usage means that the transcript occupies a larger fraction of
its gene's transcript pool after UV. Negative values indicate a smaller
fraction.

## Permutation test

Each comparison contains three control and three UV replicates. The six samples
can be divided into two groups of three in 20 ways. The script calculates the
observed difference and compares it with all 20 possible label assignments.

This gives a simple test with no distributional assumption, but the p-value
resolution is limited. In particular, the small number of permutations makes
BH correction very conservative. The raw p-value cutoff is therefore used only
to prepare a candidate list, not to claim transcriptome-wide significance.

## Output

The main tables are:

```text
results/statistical_DTE_candidate_transcripts.csv
results/statistical_DTE_candidate_protein_coding_transcripts.csv
results/statistical_DTU_candidate_transcripts.csv
results/statistical_DTU_candidate_protein_coding_transcripts.csv
```

The protein-coding subsets require a protein-coding gene and transcript type,
together with an annotated GENCODE protein ID.

## Limitations

The analysis uses TPM estimates and three replicates per condition. Isoforms
with similar sequence can be difficult to distinguish using short-read data,
and the permutation test does not model count dispersion. I therefore treat
these results as a first-pass ranking. A count-based DTU method and targeted
validation are needed for the strongest candidates.
