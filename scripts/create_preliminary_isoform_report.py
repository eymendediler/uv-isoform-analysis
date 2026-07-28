#!/usr/bin/env python3
"""
Create preliminary data tables and simple SVG figures for the isoform project.

The outputs are meant for lab discussion, thesis proposal slides, or an early
progress report. They summarize DTE and DTU candidate isoforms produced by
scripts/statistical_isoform_analysis.py.
"""

from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"
OUT_DIR = PROJECT_DIR / "preliminary_data"

DTE_PATH = RESULTS_DIR / "statistical_DTE_candidate_protein_coding_transcripts.csv"
DTU_PATH = RESULTS_DIR / "statistical_DTU_candidate_protein_coding_transcripts.csv"


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>',
        'text { font-family: Arial, Helvetica, sans-serif; fill: #222; }',
        '.axis { stroke: #555; stroke-width: 1; }',
        '.grid { stroke: #dddddd; stroke-width: 1; }',
        '.up { fill: #b83a4b; }',
        '.down { fill: #2f6f9f; }',
        '.usageup { fill: #b8792e; }',
        '.usagedown { fill: #4f7f60; }',
        '</style>',
    ]


def save_svg(path: Path, parts: list[str]) -> None:
    path.write_text("\n".join(parts + ["</svg>\n"]))


def count_table(df: pd.DataFrame, direction_col: str) -> pd.DataFrame:
    return (
        df.groupby(["time_min", direction_col])
        .size()
        .reset_index(name="count")
        .sort_values(["time_min", direction_col])
    )


def bar_chart_counts(summary: pd.DataFrame, path: Path, title: str, direction_col: str) -> None:
    width, height = 760, 430
    margin_left, margin_right, margin_top, margin_bottom = 82, 30, 62, 72
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    max_count = int(summary["count"].max())
    y_max = ((max_count // 250) + 1) * 250

    parts = svg_header(width, height)
    parts.append(f'<text x="{width/2}" y="28" text-anchor="middle" font-size="18" font-weight="700">{escape(title)}</text>')
    parts.append(f'<line x1="{margin_left}" y1="{margin_top + plot_h}" x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" class="axis"/>')
    parts.append(f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_h}" class="axis"/>')

    for tick in range(0, y_max + 1, max(250, y_max // 5)):
        y = margin_top + plot_h - (tick / y_max) * plot_h
        parts.append(f'<line x1="{margin_left}" y1="{y:.1f}" x2="{margin_left + plot_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(f'<text x="{margin_left - 10}" y="{y + 4:.1f}" text-anchor="end" font-size="12">{tick}</text>')

    times = [12, 30, 60]
    directions = list(summary[direction_col].drop_duplicates())
    group_w = plot_w / len(times)
    bar_w = 54

    for i, minute in enumerate(times):
        x_center = margin_left + group_w * (i + 0.5)
        parts.append(f'<text x="{x_center}" y="{height - 28}" text-anchor="middle" font-size="13">{minute} min</text>')
        for j, direction in enumerate(directions):
            row = summary[(summary["time_min"] == minute) & (summary[direction_col] == direction)]
            if row.empty:
                continue
            count = int(row["count"].iloc[0])
            x = x_center + (j - 0.5) * (bar_w + 16)
            bar_h = (count / y_max) * plot_h
            y = margin_top + plot_h - bar_h
            cls = "up" if direction == "up" else "down"
            if direction == "usage_up":
                cls = "usageup"
            if direction == "usage_down":
                cls = "usagedown"
            parts.append(f'<rect x="{x - bar_w/2:.1f}" y="{y:.1f}" width="{bar_w}" height="{bar_h:.1f}" class="{cls}"/>')
            parts.append(f'<text x="{x:.1f}" y="{y - 8:.1f}" text-anchor="middle" font-size="12">{count}</text>')
            parts.append(f'<text x="{x:.1f}" y="{height - 48}" text-anchor="middle" font-size="11">{escape(direction)}</text>')

    parts.append(f'<text x="24" y="{margin_top + plot_h/2}" transform="rotate(-90 24 {margin_top + plot_h/2})" text-anchor="middle" font-size="13">candidate count</text>')
    save_svg(path, parts)


def top_table(df: pd.DataFrame, value_col: str, n: int = 25) -> pd.DataFrame:
    high = df.sort_values(value_col, ascending=False).head(n).copy()
    low = df.sort_values(value_col, ascending=True).head(n).copy()
    return pd.concat([high, low], ignore_index=True)


def make_usage_examples(dtu: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "gene_name",
        "transcript_id",
        "transcript_name",
        "protein_id",
        "time_min",
        "control_mean_usage",
        "uv_mean_usage",
        "delta_usage",
        "control_gene_mean_TPM",
        "uv_gene_mean_TPM",
        "permutation_pvalue",
        "padj_BH",
    ]
    return dtu.sort_values("delta_usage", key=lambda s: s.abs(), ascending=False)[cols].head(40)


def to_markdown_table(df: pd.DataFrame) -> str:
    """Render a small DataFrame as a Markdown table without optional packages."""
    text_df = df.copy()
    for col in text_df.columns:
        if pd.api.types.is_float_dtype(text_df[col]):
            text_df[col] = text_df[col].map(lambda value: f"{value:.4g}")
        else:
            text_df[col] = text_df[col].astype(str)

    headers = list(text_df.columns)
    rows = text_df.values.tolist()
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def usage_slope_plot(examples: pd.DataFrame, path: Path) -> None:
    plot_df = examples.head(12).copy()
    width, height = 900, 560
    axis_x0, axis_x1 = 330, 810
    top, row_h = 78, 34
    parts = svg_header(width, height)
    parts.append('<text x="450" y="28" text-anchor="middle" font-size="18" font-weight="700">Top DTU protein-coding isoform candidates</text>')
    parts.append('<text x="570" y="52" text-anchor="middle" font-size="12">isoform usage fraction inside its gene</text>')
    parts.append(f'<line x1="{axis_x0}" y1="64" x2="{axis_x1}" y2="64" class="axis"/>')
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        x = axis_x0 + tick * (axis_x1 - axis_x0)
        parts.append(f'<line x1="{x:.1f}" y1="59" x2="{x:.1f}" y2="69" class="axis"/>')
        parts.append(f'<text x="{x:.1f}" y="84" text-anchor="middle" font-size="11">{tick:.2g}</text>')

    for i, row in enumerate(plot_df.itertuples(index=False)):
        y = top + i * row_h + 28
        c_x = axis_x0 + row.control_mean_usage * (axis_x1 - axis_x0)
        u_x = axis_x0 + row.uv_mean_usage * (axis_x1 - axis_x0)
        cls = "usageup" if row.delta_usage > 0 else "usagedown"
        label = f"{row.gene_name} {row.transcript_name} ({int(row.time_min)}m)"
        parts.append(f'<text x="18" y="{y + 4:.1f}" font-size="12">{escape(label[:42])}</text>')
        parts.append(f'<line x1="{axis_x0}" y1="{y:.1f}" x2="{axis_x1}" y2="{y:.1f}" class="grid"/>')
        parts.append(f'<line x1="{c_x:.1f}" y1="{y:.1f}" x2="{u_x:.1f}" y2="{y:.1f}" stroke="#777" stroke-width="1.8"/>')
        parts.append(f'<circle cx="{c_x:.1f}" cy="{y:.1f}" r="5" fill="#2f6f9f"/>')
        parts.append(f'<circle cx="{u_x:.1f}" cy="{y:.1f}" r="5" class="{cls}"/>')
        parts.append(f'<text x="{axis_x1 + 16}" y="{y + 4:.1f}" font-size="11">Δ {row.delta_usage:+.2f}</text>')
    parts.append(f'<circle cx="{axis_x0}" cy="{height - 34}" r="5" fill="#2f6f9f"/>')
    parts.append(f'<text x="{axis_x0 + 12}" y="{height - 30}" font-size="12">control</text>')
    parts.append(f'<circle cx="{axis_x0 + 90}" cy="{height - 34}" r="5" class="usageup"/>')
    parts.append(f'<text x="{axis_x0 + 102}" y="{height - 30}" font-size="12">UV usage up</text>')
    parts.append(f'<circle cx="{axis_x0 + 210}" cy="{height - 34}" r="5" class="usagedown"/>')
    parts.append(f'<text x="{axis_x0 + 222}" y="{height - 30}" font-size="12">UV usage down</text>')
    save_svg(path, parts)


def write_markdown_report(dte: pd.DataFrame, dtu: pd.DataFrame, dte_summary: pd.DataFrame, dtu_summary: pd.DataFrame) -> None:
    dte_top = dte.sort_values("log2FC", ascending=False).head(8)
    dtu_top = dtu.sort_values("delta_usage", key=lambda s: s.abs(), ascending=False).head(8)

    lines = [
        "# Preliminary Isoform Data",
        "",
        "## What We Can Say Now",
        "",
        "- We built an isoform-level RNA-seq analysis from GEO Salmon `quant.sf.gz` files.",
        "- We mapped transcript isoforms to GENCODE v35 gene names and protein IDs.",
        "- We separated two questions: DTE for expression change and DTU for isoform usage change.",
        "- These are preliminary candidates because the current exact permutation test has limited resolution with 3 vs 3 replicates.",
        "",
        "## DTE Protein-Coding Candidate Counts",
        "",
        to_markdown_table(dte_summary),
        "",
        "## DTU Protein-Coding Candidate Counts",
        "",
        to_markdown_table(dtu_summary),
        "",
        "## Strong DTE Examples",
        "",
        to_markdown_table(dte_top[["gene_name", "transcript_name", "protein_id", "time_min", "direction", "control_mean_TPM", "uv_mean_TPM", "log2FC"]]),
        "",
        "## Strong DTU Examples",
        "",
        to_markdown_table(dtu_top[["gene_name", "transcript_name", "protein_id", "time_min", "direction", "control_mean_usage", "uv_mean_usage", "delta_usage"]]),
        "",
        "## Recommended Next Steps",
        "",
        "1. Validate DTE with a count-based RNA-seq method, ideally tximport plus DESeq2 or edgeR.",
        "2. Validate DTU with a specialist method such as DRIMSeq, DEXSeq, satuRn, or IsoformSwitchAnalyzeR.",
        "3. Focus first on UV response genes and DNA damage response genes: JUN, FOS/FOSB, EGR1/EGR2, ATF3, DUSP family, NR4A1, TP53-related pathways.",
        "4. Make gene-specific isoform plots for 10-20 biologically interesting candidates.",
        "5. Check whether isoform changes are coding, NMD, retained intron, or non-coding transcript changes.",
    ]
    (OUT_DIR / "preliminary_isoform_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    dte = pd.read_csv(DTE_PATH)
    dtu = pd.read_csv(DTU_PATH)

    dte_summary = count_table(dte, "direction")
    dtu_summary = count_table(dtu, "direction")

    dte_summary.to_csv(OUT_DIR / "dte_protein_coding_candidate_counts.csv", index=False)
    dtu_summary.to_csv(OUT_DIR / "dtu_protein_coding_candidate_counts.csv", index=False)

    top_table(dte, "log2FC").to_csv(OUT_DIR / "top_dte_protein_coding_candidates.csv", index=False)
    make_usage_examples(dtu).to_csv(OUT_DIR / "top_dtu_protein_coding_candidates.csv", index=False)

    bar_chart_counts(
        dte_summary,
        OUT_DIR / "figure_1_dte_candidate_counts.svg",
        "Protein-coding DTE isoform candidates",
        "direction",
    )
    bar_chart_counts(
        dtu_summary,
        OUT_DIR / "figure_2_dtu_candidate_counts.svg",
        "Protein-coding DTU isoform candidates",
        "direction",
    )
    usage_slope_plot(
        make_usage_examples(dtu),
        OUT_DIR / "figure_3_top_dtu_usage_changes.svg",
    )

    write_markdown_report(dte, dtu, dte_summary, dtu_summary)

    print(f"Wrote preliminary data package to {OUT_DIR}")


if __name__ == "__main__":
    main()
