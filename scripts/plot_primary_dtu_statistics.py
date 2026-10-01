#!/usr/bin/env python3
"""Create replicate-aware figures for primary protein-changing DTU candidates."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results" / "dtu_v2"
OUT = RESULTS / "primary_candidate_timepoint_statistics"
TIMES = (0, 12, 30, 60)


def q_label(value: float) -> str:
    if pd.isna(value):
        return "q=NA"
    if value < 0.001:
        return "q<0.001"
    return f"q={value:.3f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(PROJECT / "metadata" / "samples.tsv", sep="\t")
    robust = pd.read_csv(
        RESULTS / "robust_protein_sequence_details_scaled_tpm_prop_sd_0_1.tsv.gz",
        sep="\t",
    )
    gene_summary = pd.read_csv(
        RESULTS / "robust_protein_sequence_gene_summary_scaled_tpm_prop_sd_0_1.tsv",
        sep="\t",
    )
    primary_genes = gene_summary.loc[
        gene_summary.multiple_distinct_protein_sequences.astype(str).str.lower()
        == "true",
        ["gene_id", "gene_name"],
    ]
    robust = robust.merge(primary_genes, on=["gene_id", "gene_name"])
    wanted_tx = set(robust.feature_id)
    wanted_genes = set(robust.gene_id)

    annotation = pd.read_csv(
        PROJECT / "metadata" / "gencode_v35_transcripts.tsv.gz", sep="\t"
    )[["transcript_id", "gene_id", "transcript_name"]]
    gene_annotation = annotation.loc[annotation.gene_id.isin(wanted_genes)]

    records = []
    for sample in samples.itertuples(index=False):
        quant = pd.read_csv(PROJECT / sample.quant_file, sep="\t", compression="gzip")
        quant = quant.merge(gene_annotation, left_on="Name", right_on="transcript_id")
        quant["gene_tpm"] = quant.groupby("gene_id").TPM.transform("sum")
        quant["usage"] = np.where(
            quant.gene_tpm > 0, quant.TPM / quant.gene_tpm, np.nan
        )
        quant = quant.loc[quant.transcript_id.isin(wanted_tx)].copy()
        quant["sample_id"] = sample.sample_id
        quant["time_min"] = sample.time_min
        quant["replicate"] = sample.replicate
        records.append(
            quant[
                [
                    "gene_id",
                    "transcript_id",
                    "transcript_name",
                    "TPM",
                    "gene_tpm",
                    "usage",
                    "sample_id",
                    "time_min",
                    "replicate",
                ]
            ]
        )
    values = pd.concat(records, ignore_index=True).merge(primary_genes, on="gene_id")

    contrasts = pd.read_csv(
        RESULTS / "timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1.tsv.gz",
        sep="\t",
    )
    contrasts = contrasts.loc[contrasts.feature_id.isin(wanted_tx)].copy()
    contrast_export = contrasts[
        [
            "gene_name",
            "feature_id",
            "transcript_name",
            "time_min",
            "tpm_mean_usage_0",
            "tpm_mean_usage_time",
            "tpm_delta_usage_vs_0",
            "direction",
            "pvalue",
            "global_bh_qvalue",
            "proportion_sd_all_samples",
        ]
    ].sort_values(["gene_name", "feature_id", "time_min"])
    contrast_export.to_csv(
        OUT / "primary_candidate_timepoint_statistics.tsv", sep="\t", index=False
    )
    values.to_csv(OUT / "primary_candidate_replicate_values.tsv", sep="\t", index=False)

    colors = plt.get_cmap("tab10")
    for gene_name, gene_values in values.groupby("gene_name", sort=True):
        gene_robust = robust.loc[robust.gene_name == gene_name].copy()
        gene_contrasts = contrasts.loc[contrasts.gene_name == gene_name]
        fig, axes = plt.subplots(
            1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [0.8, 1.6]}
        )

        # Gene abundance: the three points are the biological replicates.
        gene_abundance = gene_values.drop_duplicates(["sample_id", "gene_id"])
        for time in TIMES:
            subset = gene_abundance.loc[gene_abundance.time_min == time]
            jitter = np.linspace(-1.2, 1.2, len(subset))
            axes[0].scatter(
                time + jitter, subset.gene_tpm, color="#222222", s=34, alpha=0.8
            )
        abundance_mean = gene_abundance.groupby("time_min").gene_tpm.mean()
        axes[0].plot(
            abundance_mean.index,
            abundance_mean.values,
            color="#3366aa",
            marker="o",
            linewidth=2,
            label="mean (display only)",
        )
        axes[0].set_title("Total gene abundance")
        axes[0].set_ylabel("Gene TPM")
        axes[0].set_xlabel("Minutes after UV")
        axes[0].set_xticks(TIMES)
        axes[0].legend(frameon=False, fontsize=8)

        for tx_index, row in enumerate(gene_robust.itertuples(index=False)):
            tx_values = gene_values.loc[gene_values.transcript_id == row.feature_id]
            color = colors(tx_index % 10)
            for time in TIMES:
                subset = tx_values.loc[tx_values.time_min == time]
                jitter = np.linspace(-1.0, 1.0, len(subset))
                axes[1].scatter(
                    time + jitter,
                    subset.usage,
                    color=color,
                    s=28,
                    alpha=0.6,
                )
            means = tx_values.groupby("time_min").usage.mean().reindex(TIMES)
            label = f"{row.transcript_name} ({row.feature_id.split('.')[0]})"
            axes[1].plot(
                TIMES,
                means.values,
                color=color,
                marker="o",
                linewidth=2,
                label=label,
            )
            tx_contrasts = gene_contrasts.loc[
                gene_contrasts.feature_id == row.feature_id
            ]
            for result in tx_contrasts.itertuples(index=False):
                if result.global_bh_qvalue <= 0.05:
                    y = tx_values.loc[tx_values.time_min == result.time_min, "usage"].max()
                    axes[1].text(
                        result.time_min,
                        min(1.02, y + 0.035),
                        "*",
                        color=color,
                        ha="center",
                        va="bottom",
                        fontsize=14,
                        fontweight="bold",
                    )

        omnibus_q = gene_robust.stage_gene_padj.min()
        axes[1].set_title(
            f"Within-gene transcript usage\nomnibus stageR gene q={omnibus_q:.3g}"
        )
        axes[1].set_ylabel("Transcript TPM / gene TPM")
        axes[1].set_xlabel("Minutes after UV")
        axes[1].set_xticks(TIMES)
        axes[1].set_ylim(-0.03, 1.08)
        axes[1].legend(frameon=False, fontsize=7, loc="best")
        axes[1].text(
            0.01,
            -0.18,
            "Points = biological samples; lines = display means; * = secondary 0-min contrast, global Benjamini–Hochberg q≤0.05",
            transform=axes[1].transAxes,
            fontsize=8,
        )
        fig.suptitle(f"{gene_name}: replicate-aware DTU evidence", fontsize=14)
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        fig.savefig(OUT / f"{gene_name}_dtu_statistics.png", dpi=220)
        plt.close(fig)

    # Summary heatmap: colour is the descriptive change; stars are formal
    # secondary contrast results after one global BH correction family.
    heat = contrast_export.copy()
    heat["row_label"] = heat.gene_name + " | " + heat.transcript_name.fillna(
        heat.feature_id.str.split(".").str[0]
    )
    delta = heat.pivot(index="row_label", columns="time_min", values="tpm_delta_usage_vs_0")
    qvals = heat.pivot(index="row_label", columns="time_min", values="global_bh_qvalue")
    order = (
        heat.groupby("row_label").tpm_delta_usage_vs_0.apply(lambda x: x.abs().max())
        .sort_values(ascending=False)
        .index
    )
    delta = delta.reindex(order)[[12, 30, 60]]
    qvals = qvals.reindex(order)[[12, 30, 60]]
    lim = max(0.1, np.nanmax(np.abs(delta.to_numpy())))
    fig_h = max(7, 0.36 * len(delta) + 2.2)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    image = ax.imshow(delta.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(3), ["12 vs 0", "30 vs 0", "60 vs 0"])
    ax.set_yticks(range(len(delta)), delta.index, fontsize=8)
    for i in range(len(delta)):
        for j in range(3):
            q = qvals.iloc[i, j]
            text = f"{delta.iloc[i, j]:+.2f}"
            if pd.notna(q) and q <= 0.05:
                text += "*"
            ax.text(j, i, text, ha="center", va="center", fontsize=7)
    cbar = fig.colorbar(image, ax=ax, shrink=0.75)
    cbar.set_label("Change in within-gene TPM usage")
    ax.set_title(
        "Primary protein-changing candidates: time-localized usage shifts\n"
        "Cell value = mean usage change\n"
        "* global Benjamini–Hochberg q≤0.05 across all transcript × time tests"
    )
    fig.tight_layout()
    fig.savefig(OUT / "primary_candidate_timepoint_heatmap.png", dpi=240)
    fig.savefig(OUT / "primary_candidate_timepoint_heatmap.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
