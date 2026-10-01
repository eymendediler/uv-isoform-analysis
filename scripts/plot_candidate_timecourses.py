#!/usr/bin/env python3
"""Plot replicate-level gene abundance and transcript usage for DTU candidates."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results" / "dtu_v2"
OUT = RESULTS / "candidate_timecourses_scaled_tpm"
CANDIDATES = ("RUNX1", "ERN1", "EPC1", "NR4A1", "E2F3", "ELK4", "FOS")
TIMES = (0, 12, 30, 60)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(PROJECT / "metadata" / "samples.tsv", sep="\t")
    robust = pd.read_csv(
        RESULTS / "robust_protein_sequence_details_scaled_tpm.tsv.gz", sep="\t"
    )
    robust = robust.loc[robust.gene_name.isin(CANDIDATES)].copy()
    wanted_genes = set(robust.gene_id)

    annotation = pd.read_csv(
        PROJECT / "metadata" / "gencode_v35_transcripts.tsv.gz", sep="\t"
    )[["transcript_id", "gene_id"]]
    annotation = annotation.loc[annotation.gene_id.isin(wanted_genes)]

    records = []
    for sample in samples.itertuples(index=False):
        quant = pd.read_csv(PROJECT / sample.quant_file, sep="\t", compression="gzip")
        quant = quant.merge(annotation, left_on="Name", right_on="transcript_id")
        totals = quant.groupby("gene_id").TPM.transform("sum")
        quant["usage"] = np.where(totals > 0, quant.TPM / totals, np.nan)
        quant["sample_id"] = sample.sample_id
        quant["time_min"] = sample.time_min
        quant["replicate"] = sample.replicate
        records.append(
            quant[["gene_id", "transcript_id", "sample_id", "time_min", "replicate", "TPM", "usage"]]
        )
    long = pd.concat(records, ignore_index=True)

    gene_totals = (
        long.groupby(["gene_id", "sample_id", "time_min", "replicate"], as_index=False)
        .TPM.sum()
        .rename(columns={"TPM": "gene_total_tpm"})
    )
    selected = long.merge(
        robust[
            [
                "gene_id", "gene_name", "feature_id", "transcript_name",
                "protein_id", "protein_length_aa", "stage_tx_padj",
                "tpm_max_abs_delta_usage",
            ]
        ],
        left_on=["gene_id", "transcript_id"],
        right_on=["gene_id", "feature_id"],
        validate="many_to_one",
    )
    selected = selected.merge(
        gene_totals, on=["gene_id", "sample_id", "time_min", "replicate"]
    )
    selected.to_csv(OUT / "candidate_replicate_values.tsv", sep="\t", index=False)

    audit = (
        selected.groupby(
            ["gene_id", "gene_name", "transcript_id", "transcript_name", "protein_id"],
            as_index=False,
            dropna=False,
        )
        .agg(
            min_tpm=("TPM", "min"),
            median_tpm=("TPM", "median"),
            max_tpm=("TPM", "max"),
            max_usage=("usage", "max"),
            max_abs_delta_usage=("tpm_max_abs_delta_usage", "first"),
            stageR_adjusted_p=("stage_tx_padj", "first"),
            protein_length_aa=("protein_length_aa", "first"),
        )
    )
    audit.to_csv(OUT / "candidate_audit_summary.tsv", sep="\t", index=False)

    offsets = (-1.7, 0, 1.7)
    for gene_name in CANDIDATES:
        g = selected.loc[selected.gene_name == gene_name].copy()
        if g.empty:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
        gt = gene_totals.loc[gene_totals.gene_id == g.gene_id.iloc[0]]
        for rep, offset in zip((1, 2, 3), offsets):
            z = gt.loc[gt.replicate == rep].sort_values("time_min")
            axes[0].scatter(z.time_min + offset, z.gene_total_tpm, s=35, label=f"repeat {rep}")
        means = gt.groupby("time_min").gene_total_tpm.mean().reindex(TIMES)
        axes[0].plot(TIMES, means, color="black", linewidth=2, marker="o", label="mean")
        axes[0].set_title(f"{gene_name}: total gene abundance")
        axes[0].set_ylabel("Total gene TPM")
        axes[0].set_xlabel("Minutes after UV")
        axes[0].legend(fontsize=8)

        for tx, z in g.groupby("transcript_id"):
            name = z.transcript_name.iloc[0]
            plen = int(z.protein_length_aa.iloc[0])
            label = f"{name} ({tx.split('.')[0]}, {plen} aa)"
            for rep, offset in zip((1, 2, 3), offsets):
                q = z.loc[z.replicate == rep]
                axes[1].scatter(q.time_min + offset, q.usage, s=25, alpha=0.7)
            means = z.groupby("time_min").usage.mean().reindex(TIMES)
            axes[1].plot(TIMES, means, linewidth=2, marker="o", label=label)
        axes[1].set_title(f"{gene_name}: robust coding-isoform usage")
        axes[1].set_ylabel("Fraction of total gene TPM")
        axes[1].set_xlabel("Minutes after UV")
        axes[1].set_ylim(bottom=0)
        axes[1].legend(fontsize=7)
        for ax in axes:
            ax.set_xticks(TIMES)
            ax.grid(alpha=0.2)
        fig.suptitle("Points: biological replicates; lines: time-point means", fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / f"{gene_name}_timecourse.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
