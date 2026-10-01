#!/usr/bin/env python3
"""Summarize agreement across the three unblocked DRIMSeq/stageR filters."""

import argparse
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "results" / "dtu_v2"
FILTERS = ("loose", "medium", "strict")
ALPHA = 0.05


def load_run(filter_name: str, input_scale: str) -> pd.DataFrame:
    suffix = "_scaled_tpm" if input_scale == "scaled_tpm" else ""
    folder = ROOT / f"{filter_name}_unblocked{suffix}"
    tx = pd.read_csv(folder / "transcript_omnibus_and_usage.tsv.gz", sep="\t")
    stage_name = globals().get("STAGE_FILE", "stageR_adjusted.tsv")
    stage = pd.read_csv(folder / stage_name, sep="\t").rename(
        columns={"gene": "stage_gene_padj", "transcript": "stage_tx_padj"}
    )
    return tx.merge(
        stage,
        left_on=["gene_id", "feature_id"],
        right_on=["geneID", "txID"],
        how="left",
        validate="one_to_one",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-scale", choices=("scaled_tpm", "raw_counts"), default="scaled_tpm"
    )
    parser.add_argument("--prop-sd-filter", action="store_true")
    args = parser.parse_args()
    global STAGE_FILE
    STAGE_FILE = (
        "stageR_adjusted_prop_sd_0_1.tsv"
        if args.prop_sd_filter else "stageR_adjusted.tsv"
    )
    suffix = "_scaled_tpm" if args.input_scale == "scaled_tpm" else ""
    if args.prop_sd_filter:
        suffix += "_prop_sd_0_1"
    runs = {name: load_run(name, args.input_scale) for name in FILTERS}
    rows = []
    confirmed = {}
    for name, data in runs.items():
        confirmed[name] = set(data.loc[data.stage_tx_padj <= ALPHA, "feature_id"])
        pc = (
            (data.gene_type == "protein_coding")
            & (data.transcript_type == "protein_coding")
            & data.protein_id.notna()
        )
        rows.append(
            {
                "filter": name,
                "tested_genes": data.gene_id.nunique(),
                "tested_transcripts": len(data),
                "stageR_confirmed_genes": data.loc[
                    data.stage_tx_padj <= ALPHA, "gene_id"
                ].nunique(),
                "stageR_confirmed_transcripts": int(
                    (data.stage_tx_padj <= ALPHA).sum()
                ),
                "stageR_confirmed_protein_coding_transcripts": int(
                    ((data.stage_tx_padj <= ALPHA) & pc).sum()
                ),
            }
        )
    pd.DataFrame(rows).to_csv(
        ROOT / f"sensitivity_summary{suffix}.tsv", sep="\t", index=False
    )

    robust_ids = set.intersection(*confirmed.values())
    # Every robust transcript is retained by the strict run. TPM usage is
    # calculated against all annotated transcripts, so it is filter invariant.
    robust = runs["strict"].loc[
        runs["strict"].feature_id.isin(robust_ids)
    ].copy()
    robust["protein_coding_transcript"] = (
        (robust.gene_type == "protein_coding")
        & (robust.transcript_type == "protein_coding")
        & robust.protein_id.notna()
    )
    for cutoff in (0.05, 0.10, 0.20):
        label = str(cutoff).replace(".", "_")
        robust[f"abs_tpm_delta_ge_{label}"] = (
            robust.tpm_max_abs_delta_usage >= cutoff
        )
    robust.to_csv(
        ROOT / f"robust_stageR_transcripts_all_filters{suffix}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    pc = robust.loc[robust.protein_coding_transcript].copy()
    gene_summary = (
        pc.groupby(["gene_id", "gene_name"], dropna=False)
        .agg(
            robust_protein_coding_transcripts=("feature_id", "nunique"),
            distinct_protein_ids=("protein_id", "nunique"),
            max_abs_tpm_delta_usage=("tpm_max_abs_delta_usage", "max"),
        )
        .reset_index()
    )
    gene_summary["multiple_distinct_protein_products"] = (
        (gene_summary.robust_protein_coding_transcripts >= 2)
        & (gene_summary.distinct_protein_ids >= 2)
    )
    gene_summary.sort_values(
        ["multiple_distinct_protein_products", "max_abs_tpm_delta_usage"],
        ascending=[False, False],
    ).to_csv(
        ROOT / f"robust_protein_isoform_gene_summary{suffix}.tsv",
        sep="\t",
        index=False,
    )

    print(pd.DataFrame(rows).to_string(index=False))
    print(f"Robust stageR transcripts across all filters: {len(robust):,}")
    print(
        "Robust protein-coding transcripts: "
        f"{robust.protein_coding_transcript.sum():,}"
    )
    print(
        "Genes with >=2 robust coding transcripts and distinct protein IDs: "
        f"{gene_summary.multiple_distinct_protein_products.sum():,}"
    )


if __name__ == "__main__":
    main()
