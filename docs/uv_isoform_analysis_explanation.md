# Analysis background

## Dataset

The RNA-seq data come from Kaya and Adebali (2025), *UV-induced
reorganization of 3D genome mediates DNA damage response*. HeLa-S3 cells were
collected without UV treatment and at 12, 30 and 60 min after UV exposure.
Each condition has three biological replicates.

The GEO files contain Salmon transcript quantifications:

```text
Name             transcript ID
Length           transcript length
EffectiveLength  effective length used by Salmon
TPM              normalized transcript abundance
NumReads         estimated assigned reads
```

GENCODE v35 was used to map transcript IDs to genes, transcript biotypes and
protein IDs.

This is an RNA-level analysis. A protein ID attached to a transcript indicates
an annotated coding product, but the RNA-seq data do not show whether that
protein is present or changes in abundance.

## Initial transcript screen

The script `scripts/analyze_transcript_isoforms.py` compares each UV time point
with the shared control. For each transcript it calculates condition means,
replicate-level fold changes and:

```text
mean log2FC = log2((mean UV TPM + 0.1) /
                   (mean control TPM + 0.1))
```

The pseudocount avoids division by zero. The initial screen requires:

```text
maximum mean TPM >= 1
absolute mean log2FC >= 1
all three replicate pairs change in the same direction
```

This screen was intended to find clear, reproducible expression changes. It
does not provide a formal count-based differential expression model.

## Annotation and filtering

Results are reported both for all annotated transcripts and for a
protein-coding subset. The latter keeps transcripts for which:

```text
gene_type = protein_coding
transcript_type = protein_coding
protein_id is present
```

Versioned Ensembl transcript IDs are retained because they match the Salmon
index and GENCODE v35 annotation used in the original processing.

## Interpretation

Transcript expression and isoform usage are different measurements. A
transcript may increase because the whole gene increases, while a DTU change
means that its relative contribution within the gene has changed. For this
reason, the initial expression screen was followed by separate DTE, DTU and
SUPPA2 analyses.

The result tables are candidate lists. Three replicates provide limited power,
TPM estimates carry transcript-assignment uncertainty, and RNA-seq alone
cannot establish a change in protein isoform abundance.
