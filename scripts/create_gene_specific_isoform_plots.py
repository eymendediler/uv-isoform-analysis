#!/usr/bin/env python3
"""
Create gene-specific isoform usage plots for preliminary inspection.

Each SVG summarizes one gene across 0, 12, 30, and 60 min:
- stacked bars: mean isoform usage fraction inside the gene
- black line: total gene TPM
- small table: mean TPM values for plotted isoforms
"""

from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

from analyze_transcript_isoforms import build_expression_table, load_transcript_annotation


PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "preliminary_data" / "gene_specific_isoform_plots"

GENES = [
    "JUN",
    "EGR1",
    "EGR2",
    "NR4A1",
    "HDAC11",
    "ZNF831",
    "SPAG1",
    "SKIL",
    "RUNX1",
    "CXCL8",
]

TIMEPOINTS = [0, 12, 30, 60]
TIME_LABELS = {0: "control", 12: "12m UV", 30: "30m UV", 60: "60m UV"}
MAX_ISOFORMS_TO_PLOT = 6
MIN_MAX_MEAN_TPM = 0.25
COLORS = ["#4f7f60", "#b8792e", "#2f6f9f", "#b83a4b", "#6f5fa8", "#7f6d4f"]


def sample_cols_for_time(minute: int) -> list[str]:
    if minute == 0:
        return [f"control_rep{i}_TPM" for i in (1, 2, 3)]
    return [f"uv_{minute}min_rep{i}_TPM" for i in (1, 2, 3)]


def prepare_gene_table(expr: pd.DataFrame, annotation: pd.DataFrame, gene_name: str) -> pd.DataFrame:
    merged = annotation.merge(expr, on="transcript_id", how="inner")
    gene = merged[merged["gene_name"] == gene_name].copy()
    if gene.empty:
        return gene

    for minute in TIMEPOINTS:
        cols = sample_cols_for_time(minute)
        gene[f"mean_TPM_{minute}"] = gene[cols].mean(axis=1)

    mean_cols = [f"mean_TPM_{minute}" for minute in TIMEPOINTS]
    gene["max_mean_TPM"] = gene[mean_cols].max(axis=1)

    for minute in TIMEPOINTS:
        total = gene[f"mean_TPM_{minute}"].sum()
        if total > 0:
            gene[f"usage_{minute}"] = gene[f"mean_TPM_{minute}"] / total
        else:
            gene[f"usage_{minute}"] = 0.0

    return gene.sort_values("max_mean_TPM", ascending=False)


def simplify_isoforms(gene: pd.DataFrame) -> pd.DataFrame:
    expressed = gene[gene["max_mean_TPM"] >= MIN_MAX_MEAN_TPM].copy()
    top = expressed.head(MAX_ISOFORMS_TO_PLOT).copy()
    other = gene.loc[~gene.index.isin(top.index)].copy()
    other_total = sum(other[f"mean_TPM_{minute}"].sum() for minute in TIMEPOINTS)
    if other_total <= 0:
        return top

    # Reserve one legend/color slot for all transcripts not shown separately.
    top = expressed.head(MAX_ISOFORMS_TO_PLOT - 1).copy()
    other = gene.loc[~gene.index.isin(top.index)].copy()
    row = {
        "transcript_id": "other_expressed_isoforms",
        "gene_name": gene["gene_name"].iloc[0],
        "transcript_name": "Other expressed isoforms",
        "transcript_type": "mixed",
        "protein_id": "",
        "max_mean_TPM": other["max_mean_TPM"].max(),
    }
    for minute in TIMEPOINTS:
        row[f"mean_TPM_{minute}"] = other[f"mean_TPM_{minute}"].sum()
        row[f"usage_{minute}"] = other[f"usage_{minute}"].sum()
    top = pd.concat([top, pd.DataFrame([row])], ignore_index=True)
    return top


def write_gene_csv(gene_name: str, gene: pd.DataFrame) -> None:
    cols = [
        "gene_name",
        "transcript_id",
        "transcript_name",
        "transcript_type",
        "protein_id",
        "max_mean_TPM",
    ]
    cols += [f"mean_TPM_{minute}" for minute in TIMEPOINTS]
    cols += [f"usage_{minute}" for minute in TIMEPOINTS]
    gene[cols].to_csv(OUT_DIR / f"{gene_name}_isoform_usage_data.csv", index=False)


def svg_header(width: int, height: int, title: str, desc: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f"<title id=\"title\">{escape(title)}</title>",
        f"<desc id=\"desc\">{escape(desc)}</desc>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        "<style>",
        "text { font-family: Arial, Helvetica, sans-serif; fill: #222; }",
        ".axis { stroke: #555; stroke-width: 1; }",
        ".grid { stroke: #dddddd; stroke-width: 1; }",
        "</style>",
    ]


def make_plot(gene_name: str, gene: pd.DataFrame) -> None:
    plot = simplify_isoforms(gene)
    write_gene_csv(gene_name, gene)

    width, height = 980, 640
    bar_x0, bar_y0 = 80, 82
    bar_w, bar_h = 300, 330
    line_x0, line_y0 = 500, 82
    line_w, line_h = 380, 220
    table_x0, table_y0 = 500, 350
    times = TIMEPOINTS

    title = f"{gene_name} isoform usage and expression across UV time course"
    desc = "Stacked bars show transcript isoform usage fractions; line shows total gene TPM."
    parts = svg_header(width, height, title, desc)
    parts.append(f'<text x="490" y="32" text-anchor="middle" font-size="20" font-weight="700">{escape(title)}</text>')
    parts.append('<text x="80" y="62" font-size="13">Isoform usage fraction inside gene</text>')
    parts.append('<text x="500" y="62" font-size="13">Total gene expression</text>')

    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        y = bar_y0 + bar_h - tick * bar_h
        parts.append(f'<line x1="{bar_x0}" y1="{y:.1f}" x2="{bar_x0 + bar_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(f'<text x="{bar_x0 - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="11">{tick:.2g}</text>')
    parts.append(f'<line x1="{bar_x0}" y1="{bar_y0}" x2="{bar_x0}" y2="{bar_y0 + bar_h}" class="axis"/>')
    parts.append(f'<line x1="{bar_x0}" y1="{bar_y0 + bar_h}" x2="{bar_x0 + bar_w}" y2="{bar_y0 + bar_h}" class="axis"/>')

    group_w = bar_w / len(times)
    for i, minute in enumerate(times):
        x = bar_x0 + i * group_w + 18
        y_current = bar_y0 + bar_h
        for j, row in enumerate(plot.itertuples(index=False)):
            usage = getattr(row, f"usage_{minute}")
            segment_h = usage * bar_h
            y_current -= segment_h
            parts.append(
                f'<rect x="{x:.1f}" y="{y_current:.1f}" width="42" height="{segment_h:.1f}" fill="{COLORS[j % len(COLORS)]}"/>'
            )
        parts.append(f'<text x="{x + 21:.1f}" y="{bar_y0 + bar_h + 22}" text-anchor="middle" font-size="12">{TIME_LABELS[minute]}</text>')

    legend_y = 450
    for j, row in enumerate(plot.itertuples(index=False)):
        y = legend_y + j * 24
        label = f"{row.transcript_name} ({row.transcript_type})"
        parts.append(f'<rect x="80" y="{y - 10}" width="12" height="12" fill="{COLORS[j % len(COLORS)]}"/>')
        parts.append(f'<text x="100" y="{y}" font-size="12">{escape(label[:48])}</text>')

    totals = [gene[f"mean_TPM_{minute}"].sum() for minute in times]
    max_total = max(totals) if max(totals) > 0 else 1
    for tick in [0, 0.5, 1.0]:
        y = line_y0 + line_h - tick * line_h
        value = tick * max_total
        parts.append(f'<line x1="{line_x0}" y1="{y:.1f}" x2="{line_x0 + line_w}" y2="{y:.1f}" class="grid"/>')
        parts.append(f'<text x="{line_x0 - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="11">{value:.1f}</text>')
    parts.append(f'<line x1="{line_x0}" y1="{line_y0}" x2="{line_x0}" y2="{line_y0 + line_h}" class="axis"/>')
    parts.append(f'<line x1="{line_x0}" y1="{line_y0 + line_h}" x2="{line_x0 + line_w}" y2="{line_y0 + line_h}" class="axis"/>')

    points = []
    for i, minute in enumerate(times):
        x = line_x0 + i * (line_w / (len(times) - 1))
        y = line_y0 + line_h - (totals[i] / max_total) * line_h
        points.append((x, y, totals[i], minute))
    path_points = " ".join([f"{x:.1f},{y:.1f}" for x, y, _, _ in points])
    parts.append(f'<polyline points="{path_points}" fill="none" stroke="#222222" stroke-width="2"/>')
    for x, y, value, minute in points:
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#222222"/>')
        parts.append(f'<text x="{x:.1f}" y="{y - 10:.1f}" text-anchor="middle" font-size="11">{value:.1f}</text>')
        parts.append(f'<text x="{x:.1f}" y="{line_y0 + line_h + 22}" text-anchor="middle" font-size="12">{TIME_LABELS[minute]}</text>')

    parts.append(f'<text x="{table_x0}" y="{table_y0}" font-size="13">Mean TPM by isoform</text>')
    header_y = table_y0 + 24
    parts.append(f'<text x="{table_x0}" y="{header_y}" font-size="11">isoform</text>')
    for i, minute in enumerate(times):
        parts.append(f'<text x="{table_x0 + 190 + i * 62}" y="{header_y}" font-size="11">{TIME_LABELS[minute]}</text>')
    for j, row in enumerate(plot.itertuples(index=False)):
        y = header_y + 22 + j * 22
        parts.append(f'<text x="{table_x0}" y="{y}" font-size="11">{escape(str(row.transcript_name)[:26])}</text>')
        for i, minute in enumerate(times):
            value = getattr(row, f"mean_TPM_{minute}")
            parts.append(f'<text x="{table_x0 + 190 + i * 62}" y="{y}" font-size="11">{value:.2f}</text>')

    output = OUT_DIR / f"{gene_name}_isoform_usage_plot.svg"
    output.write_text("\n".join(parts + ["</svg>\n"]))


def write_index(generated: list[str]) -> None:
    lines = [
        "# Gene-Specific Isoform Usage Plots",
        "",
        "Each gene plot contains:",
        "",
        "- stacked bars: isoform usage fraction inside the gene; transcripts not shown separately are combined as `Other expressed isoforms`",
        "- black line: total gene TPM across the time course",
        "- table: mean TPM of the plotted isoforms",
        "",
        "Important: these plots show mean values across 3 replicates per condition. They are for preliminary visual inspection and candidate prioritization.",
        "",
    ]
    for gene_name in generated:
        lines.append(f"- `{gene_name}_isoform_usage_plot.svg`")
        lines.append(f"- `{gene_name}_isoform_usage_data.csv`")
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expr = build_expression_table()
    annotation = load_transcript_annotation()

    generated = []
    for gene_name in GENES:
        gene = prepare_gene_table(expr, annotation, gene_name)
        if gene.empty:
            print(f"Skipping {gene_name}: no transcripts found")
            continue
        make_plot(gene_name, gene)
        generated.append(gene_name)
        print(f"Wrote {gene_name}")

    write_index(generated)
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
