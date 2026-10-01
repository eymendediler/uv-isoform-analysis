#!/usr/bin/env python3
"""Check whether robust coding transcripts encode distinct amino-acid strings."""

from __future__ import annotations

import gzip
import hashlib
import argparse
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
FASTA = PROJECT / "gencode.v35.pc_translations.fa.gz"


def read_fasta() -> dict[str, str]:
    sequences: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with gzip.open(FASTA, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks)
                current = line[1:].split("|", 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if current:
        sequences[current] = "".join(chunks)
    return sequences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-scale", choices=("scaled_tpm", "raw_counts"), default="scaled_tpm"
    )
    parser.add_argument("--prop-sd-filter", action="store_true")
    args = parser.parse_args()
    suffix = "_scaled_tpm" if args.input_scale == "scaled_tpm" else ""
    if args.prop_sd_filter:
        suffix += "_prop_sd_0_1"
    root = PROJECT / "results" / "dtu_v2"
    dtu = root / f"robust_stageR_transcripts_all_filters{suffix}.tsv.gz"
    out_detail = root / f"robust_protein_sequence_details{suffix}.tsv.gz"
    out_genes = root / f"robust_protein_sequence_gene_summary{suffix}.tsv"
    robust = pd.read_csv(dtu, sep="\t")
    coding = robust.loc[robust.protein_coding_transcript].copy()
    sequences = read_fasta()
    coding["amino_acid_sequence"] = coding.protein_id.map(sequences)
    if coding.amino_acid_sequence.isna().any():
        missing = coding.loc[coding.amino_acid_sequence.isna(), "protein_id"].tolist()
        raise RuntimeError(f"Missing GENCODE protein sequences: {missing[:5]}")
    coding["protein_length_aa"] = coding.amino_acid_sequence.str.len()
    coding["protein_sequence_sha256"] = coding.amino_acid_sequence.map(
        lambda value: hashlib.sha256(value.encode()).hexdigest()
    )
    coding.drop(columns="amino_acid_sequence").to_csv(
        out_detail, sep="\t", index=False, compression="gzip"
    )

    genes = (
        coding.groupby(["gene_id", "gene_name"], dropna=False)
        .agg(
            robust_coding_transcripts=("feature_id", "nunique"),
            distinct_protein_ids=("protein_id", "nunique"),
            distinct_amino_acid_sequences=("protein_sequence_sha256", "nunique"),
            min_protein_length_aa=("protein_length_aa", "min"),
            max_protein_length_aa=("protein_length_aa", "max"),
            max_abs_tpm_delta_usage=("tpm_max_abs_delta_usage", "max"),
        )
        .reset_index()
    )
    genes["multiple_distinct_protein_sequences"] = (
        (genes.robust_coding_transcripts >= 2)
        & (genes.distinct_amino_acid_sequences >= 2)
    )
    genes.sort_values(
        ["multiple_distinct_protein_sequences", "max_abs_tpm_delta_usage"],
        ascending=[False, False],
    ).to_csv(out_genes, sep="\t", index=False)
    print(
        "Genes with >=2 robust coding transcripts and genuinely distinct "
        f"amino-acid sequences: {genes.multiple_distinct_protein_sequences.sum():,}"
    )


if __name__ == "__main__":
    main()
