#!/usr/bin/env python3
"""Integrate robust DTU, replicate variability, protein sequence and gene DE."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
DTU = PROJECT / "results" / "dtu_v2"
CANDIDATES = ("RUNX1", "ERN1", "EPC1", "NR4A1", "E2F3", "ELK4", "FOS")
TIMES = (0, 12, 30, 60)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prop-sd-filter", action="store_true")
    args = parser.parse_args()
    suffix = "_prop_sd_0_1" if args.prop_sd_filter else ""
    out = DTU / f"biological_evidence_scaled_tpm{suffix}"
    detail_file = DTU / f"robust_protein_sequence_details_scaled_tpm{suffix}.tsv.gz"
    gene_file = DTU / f"robust_protein_sequence_gene_summary_scaled_tpm{suffix}.tsv"
    out.mkdir(parents=True, exist_ok=True)
    detail = pd.read_csv(detail_file, sep="\t")
    gene_seq = pd.read_csv(gene_file, sep="\t")
    gene_seq = gene_seq.loc[gene_seq.multiple_distinct_protein_sequences].copy()
    de = pd.read_csv(
        PROJECT / "results" / "gene_expression_time" / "gene_time_omnibus_deseq2.tsv.gz",
        sep="\t",
    )
    reps = pd.read_csv(
        DTU / "candidate_timecourses_scaled_tpm" / "candidate_replicate_values.tsv",
        sep="\t",
    )

    gene_level = gene_seq.merge(
        de[
            [
                "gene_id", "padj", "baseMean", "max_abs_descriptive_log2fc",
                "descriptive_log2fc_12_vs_0", "descriptive_log2fc_30_vs_0",
                "descriptive_log2fc_60_vs_0",
            ]
        ],
        on="gene_id",
        how="left",
        validate="one_to_one",
    )
    gene_level["gene_time_de_fdr_0_05"] = gene_level.padj <= 0.05

    med_tpm = (
        detail.groupby(["gene_id", "gene_name", "feature_id"], as_index=False)
        .agg(transcript_effect=("tpm_max_abs_delta_usage", "first"))
    )
    # Candidate replicate values provide actual transcript TPM summaries for the
    # seven hand-audited genes; leave other genes unscored rather than inventing
    # an abundance threshold from time-point means.
    cand_abundance = (
        reps.groupby(["gene_id", "gene_name", "transcript_id"], as_index=False)
        .agg(median_transcript_tpm=("TPM", "median"))
        .groupby(["gene_id", "gene_name"], as_index=False)
        .agg(lowest_median_robust_tx_tpm=("median_transcript_tpm", "min"))
    )
    del med_tpm
    gene_level = gene_level.merge(cand_abundance, on=["gene_id", "gene_name"], how="left")
    gene_level.to_csv(out / "distinct_protein_isoform_gene_evidence.tsv", sep="\t", index=False)

    cand = detail.loc[detail.gene_name.isin(CANDIDATES)].copy()
    time_sd = (
        reps.groupby(["gene_id", "gene_name", "transcript_id", "time_min"])
        .usage.std()
        .groupby(["gene_id", "gene_name"])
        .max()
        .rename("max_within_time_usage_sd")
        .reset_index()
    )
    cand_summary = (
        cand.groupby(["gene_id", "gene_name"], as_index=False)
        .agg(
            robust_coding_transcripts=("feature_id", "nunique"),
            distinct_protein_sequences=("protein_sequence_sha256", "nunique"),
            max_abs_usage_change=("tpm_max_abs_delta_usage", "max"),
            min_protein_length=("protein_length_aa", "min"),
            max_protein_length=("protein_length_aa", "max"),
        )
        .merge(cand_abundance, on=["gene_id", "gene_name"])
        .merge(time_sd, on=["gene_id", "gene_name"])
        .merge(
            de[
                ["gene_id", "padj", "baseMean", "max_abs_descriptive_log2fc",
                 "descriptive_log2fc_12_vs_0", "descriptive_log2fc_30_vs_0",
                 "descriptive_log2fc_60_vs_0"]
            ],
            on="gene_id",
        )
    )
    cand_summary["usage_effect_to_max_within_time_sd"] = (
        cand_summary.max_abs_usage_change / cand_summary.max_within_time_usage_sd
    )
    cand_summary["gene_time_de_fdr_0_05"] = cand_summary.padj <= 0.05
    cand_summary.sort_values(
        ["usage_effect_to_max_within_time_sd", "max_abs_usage_change"], ascending=False
    ).to_csv(out / "candidate_integrated_evidence.tsv", sep="\t", index=False)

    # Plot 1: separate overall gene abundance response from isoform redistribution.
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    color = -np.log10(cand_summary.padj.clip(lower=1e-300))
    sizes = 60 + 140 * np.clip(
        np.log10(cand_summary.lowest_median_robust_tx_tpm + 1), 0, 1.5
    ) / 1.5
    sc = ax.scatter(
        cand_summary.max_abs_descriptive_log2fc,
        cand_summary.max_abs_usage_change,
        c=color,
        s=sizes,
        cmap="viridis",
        edgecolor="black",
    )
    for row in cand_summary.itertuples():
        ax.annotate(row.gene_name, (row.max_abs_descriptive_log2fc, row.max_abs_usage_change),
                    xytext=(5, 4), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Maximum absolute gene-abundance change (descriptive |log2FC|)")
    ax.set_ylabel("Maximum absolute transcript-usage change")
    ax.set_title("Gene-level UV response versus coding-isoform redistribution")
    ax.grid(alpha=0.2)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Gene time-test strength (-log10 adjusted p)")
    fig.tight_layout()
    fig.savefig(out / "candidate_gene_expression_vs_isoform_change.png", dpi=200)
    plt.close(fig)

    # Plot 2: temporal usage heatmap for the most-changing robust coding
    # transcript in each distinct-protein candidate gene.
    top = (
        detail.loc[detail.gene_id.isin(gene_seq.gene_id)]
        .sort_values("tpm_max_abs_delta_usage", ascending=False)
        .drop_duplicates("gene_id")
        .copy()
    )
    cols = [f"tpm_mean_usage_{t}" for t in TIMES]
    top["peak_time"] = np.asarray(TIMES)[np.argmax(top[cols].to_numpy(), axis=1)]
    top = top.sort_values(["peak_time", "tpm_max_abs_delta_usage"], ascending=[True, False])
    matrix = top[cols].to_numpy()
    height = max(8, 0.32 * len(top))
    fig, ax = plt.subplots(figsize=(7.5, height))
    im = ax.imshow(matrix, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(4), [str(t) for t in TIMES])
    ax.set_xlabel("Minutes after UV")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(
        [f"{r.gene_name} | {r.transcript_name}" for r in top.itertuples()], fontsize=7
    )
    ax.set_title("Most-changing robust coding transcript per distinct-protein gene")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cb.set_label("Fraction of total gene TPM")
    fig.tight_layout()
    fig.savefig(out / "distinct_protein_candidate_usage_heatmap.png", dpi=220)
    plt.close(fig)

    # Plot 3: effect size versus within-time replicate spread for audited genes.
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(
        cand_summary.max_within_time_usage_sd,
        cand_summary.max_abs_usage_change,
        s=90,
        color="#3b82f6",
        edgecolor="black",
    )
    for row in cand_summary.itertuples():
        ax.annotate(row.gene_name, (row.max_within_time_usage_sd, row.max_abs_usage_change),
                    xytext=(5, 4), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Largest within-time SD across biological replicates")
    ax.set_ylabel("Maximum absolute transcript-usage change")
    ax.set_title("Observed usage effect versus replicate spread")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out / "candidate_effect_vs_replicate_spread.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
