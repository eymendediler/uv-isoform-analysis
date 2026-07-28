#!/usr/bin/env python3
"""Create publication-style overview and candidate isoform figures."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_transcript_isoforms import build_expression_table, load_transcript_annotation


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"
OUT_DIR = PROJECT_DIR / "preliminary_data" / "publication_figures"

TIMES = [0, 12, 30, 60]
TIME_LABELS = ["Control", "12 min", "30 min", "60 min"]
COLORS = {
    "down": "#0072B2",
    "up": "#D55E00",
    "overlap": "#009E73",
    "suppa": "#56B4E9",
    "permutation": "#E69F00",
}
GENE_TRANSCRIPTS = {
    "ZNF831": ["ZNF831-201", "ZNF831-203"],
    "RUNX1": ["RUNX1-202", "RUNX1-203"],
    "SKIL": ["SKIL-201", "SKIL-205"],
    "HDAC11": ["HDAC11-201", "HDAC11-213"],
}


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.8, length=3)


def label_bars(ax: plt.Axes, bars) -> None:
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{int(height):,}",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )


def create_overview_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dte = pd.read_csv(
        RESULTS_DIR / "statistical_DTE_protein_coding_candidates_summary.csv"
    )
    dte["analysis"] = "DTE"
    dte["direction"] = dte["direction"].replace({"down": "down", "up": "up"})

    dtu = pd.read_csv(
        RESULTS_DIR / "statistical_DTU_protein_coding_candidates_summary.csv"
    )
    dtu["analysis"] = "DTU"
    dtu["direction"] = dtu["direction"].replace(
        {"usage_down": "down", "usage_up": "up"}
    )

    counts = pd.concat(
        [
            dte[["analysis", "time_min", "direction", "count"]],
            dtu[["analysis", "time_min", "direction", "count"]],
        ],
        ignore_index=True,
    )
    concordance = pd.read_csv(
        RESULTS_DIR / "suppa2" / "method_concordance_summary.csv"
    )
    local = pd.read_csv(RESULTS_DIR / "suppa2" / "suppa2_candidate_summary.csv")
    local = local.loc[local["analysis"].eq("local_event")].copy()
    return counts, concordance, local


def create_overview_figure() -> None:
    counts, concordance, local = create_overview_data()
    counts.to_csv(OUT_DIR / "figure_4_overview_source_data.csv", index=False)
    concordance.to_csv(
        OUT_DIR / "figure_4_concordance_source_data.csv", index=False
    )
    local.to_csv(OUT_DIR / "figure_4_local_events_source_data.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))
    x = np.arange(3)
    width = 0.34

    for ax, analysis, panel, title in [
        (axes[0, 0], "DTE", "A", "Protein-coding transcript expression"),
        (axes[0, 1], "DTU", "B", "Protein-coding transcript usage"),
    ]:
        table = (
            counts.loc[counts["analysis"].eq(analysis)]
            .pivot(index="time_min", columns="direction", values="count")
            .reindex([12, 30, 60])
        )
        bars_down = ax.bar(
            x - width / 2,
            table["down"],
            width,
            color=COLORS["down"],
            label="Decreased",
        )
        bars_up = ax.bar(
            x + width / 2,
            table["up"],
            width,
            color=COLORS["up"],
            label="Increased",
        )
        label_bars(ax, bars_down)
        label_bars(ax, bars_up)
        ax.set_xticks(x, ["12", "30", "60"])
        ax.set_xlabel("Time after UV (min)")
        ax.set_ylabel("Candidate transcripts")
        ax.set_title(title, pad=32)
        ax.legend(
            frameon=False,
            ncol=2,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.01),
            borderaxespad=0,
            columnspacing=2.0,
            handletextpad=0.7,
        )
        ax.text(-0.14, 1.15, panel, transform=ax.transAxes, fontweight="bold", fontsize=11)
        clean_axis(ax)

    ax = axes[1, 0]
    pct_perm = 100 * concordance["overlap"] / concordance["permutation_candidates"]
    pct_suppa = 100 * concordance["overlap"] / concordance["suppa2_candidates"]
    b1 = ax.bar(
        x - width / 2,
        pct_perm,
        width,
        color=COLORS["permutation"],
        label="Of permutation candidates",
    )
    b2 = ax.bar(
        x + width / 2,
        pct_suppa,
        width,
        color=COLORS["suppa"],
        label="Of SUPPA2 candidates",
    )
    for bars in (b1, b2):
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{bar.get_height():.0f}%",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )
    for i, row in enumerate(concordance.itertuples()):
        ax.text(i, 4, f"n={row.overlap:,}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x, ["12", "30", "60"])
    ax.set_xlabel("Time after UV (min)")
    ax.set_ylabel("Cross-method overlap (%)")
    ax.set_ylim(0, 120)
    ax.set_title("Agreement between DTU methods", pad=8)
    ax.legend(frameon=False, fontsize=7, loc="upper center")
    ax.text(-0.14, 1.08, "C", transform=ax.transAxes, fontweight="bold", fontsize=11)
    clean_axis(ax)

    ax = axes[1, 1]
    event_order = ["SE", "A3", "A5", "MX", "RI", "AF", "AL"]
    heat = (
        local.pivot(index="event_type", columns="time_min", values="candidate_count")
        .reindex(event_order)
        .reindex(columns=[12, 30, 60])
    )
    image = ax.imshow(heat, cmap="Blues", aspect="auto")
    ax.set_xticks(np.arange(3), ["12", "30", "60"])
    ax.set_yticks(np.arange(len(event_order)), event_order)
    ax.set_xlabel("Time after UV (min)")
    ax.set_ylabel("Splicing event")
    ax.set_title("SUPPA2 local splicing candidates", pad=8)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            value = int(heat.iloc[i, j])
            color = "white" if value > heat.to_numpy().max() * 0.55 else "#222222"
            ax.text(j, i, f"{value:,}", ha="center", va="center", color=color, fontsize=7.5)
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Candidate events", fontsize=8)
    ax.text(-0.14, 1.08, "D", transform=ax.transAxes, fontweight="bold", fontsize=11)

    fig.suptitle(
        "UV-associated transcript and splicing changes",
        fontsize=12,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        -0.01,
        "Candidate-level results; three biological replicates per condition.",
        ha="center",
        fontsize=7.5,
    )
    fig.tight_layout(h_pad=2.4, w_pad=1.8)
    save_figure(fig, "figure_4_analysis_overview")


def sample_columns(minute: int) -> list[str]:
    prefix = "control" if minute == 0 else f"uv_{minute}min"
    return [f"{prefix}_rep{i}_TPM" for i in (1, 2, 3)]


def create_candidate_usage_data() -> pd.DataFrame:
    expr = build_expression_table().set_index("transcript_id")
    annotation = load_transcript_annotation().set_index("transcript_id")
    gene_for_transcript = annotation["gene_name"].reindex(expr.index)

    records = []
    for gene, transcript_names in GENE_TRANSCRIPTS.items():
        gene_ids = gene_for_transcript.index[gene_for_transcript.eq(gene)]
        gene_expr = expr.loc[gene_ids]
        name_map = annotation.loc[gene_ids, "transcript_name"]
        selected = name_map.index[name_map.isin(transcript_names)]
        if len(selected) != len(transcript_names):
            raise ValueError(f"Could not find all selected transcripts for {gene}")

        for minute, time_label in zip(TIMES, TIME_LABELS):
            cols = sample_columns(minute)
            totals = gene_expr[cols].sum(axis=0)
            for transcript_id in selected:
                for replicate, col in enumerate(cols, start=1):
                    usage = (
                        gene_expr.loc[transcript_id, col] / totals[col]
                        if totals[col] > 0
                        else 0.0
                    )
                    records.append(
                        {
                            "gene_name": gene,
                            "transcript_name": name_map.loc[transcript_id],
                            "transcript_id": transcript_id,
                            "time_min": minute,
                            "condition": time_label,
                            "replicate": replicate,
                            "TPM": gene_expr.loc[transcript_id, col],
                            "gene_total_TPM": totals[col],
                            "usage": usage,
                        }
                    )
    return pd.DataFrame(records)


def create_candidate_figure() -> None:
    data = create_candidate_usage_data()
    data.to_csv(OUT_DIR / "figure_5_candidate_usage_source_data.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), sharex=True, sharey=True)
    transcript_colors = ["#0072B2", "#D55E00"]
    offsets = [-0.08, 0.08]
    replicate_jitter = [-0.035, 0.0, 0.035]

    for panel_index, (ax, (gene, transcript_names)) in enumerate(
        zip(axes.flat, GENE_TRANSCRIPTS.items())
    ):
        gene_data = data.loc[data["gene_name"].eq(gene)]
        for k, transcript_name in enumerate(transcript_names):
            tx = gene_data.loc[gene_data["transcript_name"].eq(transcript_name)]
            means = tx.groupby("time_min")["usage"].mean().reindex(TIMES)
            mean_x = np.arange(4) + offsets[k]
            ax.plot(
                mean_x,
                means,
                color=transcript_colors[k],
                marker="o",
                linewidth=1.8,
                markersize=4.5,
                label=transcript_name,
                zorder=3,
            )
            for time_index, minute in enumerate(TIMES):
                values = (
                    tx.loc[tx["time_min"].eq(minute)]
                    .sort_values("replicate")["usage"]
                    .to_numpy()
                )
                ax.scatter(
                    time_index + offsets[k] + np.array(replicate_jitter),
                    values,
                    s=17,
                    facecolors="white",
                    edgecolors=transcript_colors[k],
                    linewidths=0.9,
                    zorder=4,
                )

        ax.set_title(gene, fontweight="bold")
        ax.set_xticks(np.arange(4), TIME_LABELS)
        ax.set_ylim(-0.04, 1.04)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.5)
        ax.legend(frameon=False, loc="best")
        ax.text(
            -0.13,
            1.07,
            chr(ord("A") + panel_index),
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=11,
        )
        clean_axis(ax)

    axes[0, 0].set_ylabel("Isoform usage fraction")
    axes[1, 0].set_ylabel("Isoform usage fraction")
    axes[1, 0].set_xlabel("Time after UV")
    axes[1, 1].set_xlabel("Time after UV")
    fig.suptitle(
        "Replicate-level isoform usage in selected candidate genes",
        fontsize=12,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        -0.015,
        "Open circles show individual biological replicates (n=3); filled points and lines show means.",
        ha="center",
        fontsize=7.5,
    )
    fig.tight_layout(h_pad=1.8, w_pad=1.5)
    save_figure(fig, "figure_5_candidate_isoform_usage")


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.svg")
    fig.savefig(OUT_DIR / f"{stem}.pdf")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300)
    plt.close(fig)


def validate_outputs() -> None:
    """Fail if figure source data or exported files are incomplete."""
    usage = pd.read_csv(OUT_DIR / "figure_5_candidate_usage_source_data.csv")
    group_sizes = usage.groupby(
        ["gene_name", "transcript_id", "time_min"], dropna=False
    ).size()
    if not group_sizes.eq(3).all():
        raise AssertionError("Each candidate/time point must contain three replicates")
    replicate_sets = usage.groupby(["gene_name", "transcript_id", "time_min"])[
        "replicate"
    ].agg(lambda x: set(x))
    if not replicate_sets.eq({1, 2, 3}).all():
        raise AssertionError("Replicate identifiers must be 1, 2 and 3")
    expected_usage = np.divide(
        usage["TPM"],
        usage["gene_total_TPM"],
        out=np.zeros(len(usage), dtype=float),
        where=usage["gene_total_TPM"].gt(0),
    )
    if not np.allclose(usage["usage"], expected_usage, atol=1e-12):
        raise AssertionError("Candidate usage values do not match TPM/gene TPM")
    if not usage["usage"].between(0, 1).all():
        raise AssertionError("Isoform usage must remain between 0 and 1")

    expected_files = [
        OUT_DIR / f"{stem}.{extension}"
        for stem in (
            "figure_4_analysis_overview",
            "figure_5_candidate_isoform_usage",
        )
        for extension in ("svg", "pdf", "png")
    ]
    if not all(path.exists() and path.stat().st_size > 0 for path in expected_files):
        raise AssertionError("One or more figure exports are missing or empty")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    set_style()
    create_overview_figure()
    create_candidate_figure()
    validate_outputs()
    print(f"Wrote publication figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
