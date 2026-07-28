#!/usr/bin/env python3
"""
Create validation-oriented outputs for isoform candidates.

Outputs:
- replicate dot plots for selected genes
- CDS/protein-length annotation tables for selected genes and DTU/DTE candidates

These outputs help answer:
1. Are isoform changes consistent across replicates?
2. Could the changed transcript isoforms alter the encoded protein/ORF?
"""

from pathlib import Path
from xml.sax.saxutils import escape
import gzip
import re

import numpy as np
import pandas as pd

from analyze_transcript_isoforms import (
    GENCODE_GTF,
    build_expression_table,
    load_transcript_annotation,
)


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"
OUT_DIR = PROJECT_DIR / "preliminary_data" / "isoform_validation"

GENES = ["EGR2", "HDAC11", "JUN", "EGR1", "NR4A1", "ZNF831", "SPAG1", "SKIL", "RUNX1", "CXCL8"]
TIMEPOINTS = [0, 12, 30, 60]
TIME_LABELS = {0: "control", 12: "12m", 30: "30m", 60: "60m"}
MAX_TRANSCRIPTS_PER_GENE = 4
PSEUDOCOUNT = 0.1


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', attribute_text))


def sample_cols_for_time(minute: int) -> list[str]:
    if minute == 0:
        return [f"control_rep{i}_TPM" for i in (1, 2, 3)]
    return [f"uv_{minute}min_rep{i}_TPM" for i in (1, 2, 3)]


def load_cds_annotation() -> pd.DataFrame:
    """Parse CDS lengths from GENCODE GTF."""
    rows = []
    with gzip.open(GENCODE_GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            rows.append(
                {
                    "transcript_id": attrs.get("transcript_id", ""),
                    "protein_id": attrs.get("protein_id", ""),
                    "cds_segment_bp": end - start + 1,
                }
            )

    if not rows:
        return pd.DataFrame(columns=["transcript_id", "protein_id", "cds_bp", "cds_segments", "protein_length_aa_est"])

    cds = pd.DataFrame(rows)
    out = (
        cds.groupby(["transcript_id", "protein_id"], dropna=False)
        .agg(cds_bp=("cds_segment_bp", "sum"), cds_segments=("cds_segment_bp", "size"))
        .reset_index()
    )
    out["protein_length_aa_est"] = (out["cds_bp"] / 3).round(1)
    return out


def prepare_gene_data(expr: pd.DataFrame, annotation: pd.DataFrame, gene_name: str) -> pd.DataFrame:
    gene = annotation.merge(expr, on="transcript_id", how="inner")
    gene = gene[gene["gene_name"] == gene_name].copy()
    if gene.empty:
        return gene

    for minute in TIMEPOINTS:
        cols = sample_cols_for_time(minute)
        gene[f"mean_TPM_{minute}"] = gene[cols].mean(axis=1)
    mean_cols = [f"mean_TPM_{minute}" for minute in TIMEPOINTS]
    gene["max_mean_TPM"] = gene[mean_cols].max(axis=1)
    return gene.sort_values("max_mean_TPM", ascending=False)


def svg_header(width: int, height: int, title: str, desc: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f'<title id="title">{escape(title)}</title>',
        f'<desc id="desc">{escape(desc)}</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        "<style>",
        "text { font-family: Arial, Helvetica, sans-serif; fill: #222; }",
        ".axis { stroke: #555; stroke-width: 1; }",
        ".grid { stroke: #dddddd; stroke-width: 1; }",
        ".mean { stroke: #222222; stroke-width: 1.8; fill: none; }",
        ".dot { fill: #2f6f9f; opacity: 0.78; }",
        "</style>",
    ]


def create_replicate_dot_plot(gene_name: str, gene: pd.DataFrame) -> None:
    plot = gene.head(MAX_TRANSCRIPTS_PER_GENE).copy()
    width = 980
    panel_h = 145
    height = 88 + panel_h * len(plot) + 34
    x0, x1 = 220, 900

    title = f"{gene_name} replicate-level isoform TPM"
    desc = "Dots show replicate TPM values; line connects condition means on log2 TPM scale."
    parts = svg_header(width, height, title, desc)
    parts.append(f'<text x="490" y="30" text-anchor="middle" font-size="20" font-weight="700">{escape(title)}</text>')
    parts.append('<text x="490" y="54" text-anchor="middle" font-size="12">y-axis is log2(TPM + 0.1); dots are individual replicates, line is condition mean</text>')

    x_positions = {minute: x0 + i * ((x1 - x0) / (len(TIMEPOINTS) - 1)) for i, minute in enumerate(TIMEPOINTS)}
    jitter = [-10, 0, 10]

    for idx, row in enumerate(plot.itertuples(index=False)):
        y_top = 82 + idx * panel_h
        y_bottom = y_top + 95
        transcript_label = f"{row.transcript_name} | {row.transcript_type}"
        parts.append(f'<text x="18" y="{y_top + 20}" font-size="12">{escape(transcript_label[:30])}</text>')
        parts.append(f'<text x="18" y="{y_top + 40}" font-size="11">{escape(str(row.transcript_id))}</text>')

        values = []
        means = []
        for minute in TIMEPOINTS:
            cols = sample_cols_for_time(minute)
            reps = [getattr(row, col) for col in cols]
            values.extend(reps)
            means.append(np.mean(reps))
        log_values = [np.log2(v + PSEUDOCOUNT) for v in values]
        y_min = min(log_values + [np.log2(PSEUDOCOUNT)])
        y_max = max(log_values + [np.log2(PSEUDOCOUNT)])
        if y_max == y_min:
            y_max += 1

        for tick_fraction in [0, 0.5, 1]:
            y = y_bottom - tick_fraction * (y_bottom - y_top)
            parts.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x1}" y2="{y:.1f}" class="grid"/>')
        parts.append(f'<line x1="{x0}" y1="{y_bottom}" x2="{x1}" y2="{y_bottom}" class="axis"/>')
        parts.append(f'<line x1="{x0}" y1="{y_top}" x2="{x0}" y2="{y_bottom}" class="axis"/>')

        mean_points = []
        for i, minute in enumerate(TIMEPOINTS):
            cols = sample_cols_for_time(minute)
            reps = [getattr(row, col) for col in cols]
            x = x_positions[minute]
            for rep_idx, value in enumerate(reps):
                log_value = np.log2(value + PSEUDOCOUNT)
                y = y_bottom - ((log_value - y_min) / (y_max - y_min)) * (y_bottom - y_top)
                parts.append(f'<circle cx="{x + jitter[rep_idx]:.1f}" cy="{y:.1f}" r="4" class="dot"/>')
            mean_value = np.mean(reps)
            mean_log = np.log2(mean_value + PSEUDOCOUNT)
            mean_y = y_bottom - ((mean_log - y_min) / (y_max - y_min)) * (y_bottom - y_top)
            mean_points.append((x, mean_y, mean_value))
            parts.append(f'<text x="{x}" y="{y_bottom + 20}" text-anchor="middle" font-size="11">{TIME_LABELS[minute]}</text>')
            parts.append(f'<text x="{x}" y="{mean_y - 8:.1f}" text-anchor="middle" font-size="10">{mean_value:.2f}</text>')
        path_points = " ".join(f"{x:.1f},{y:.1f}" for x, y, _ in mean_points)
        parts.append(f'<polyline points="{path_points}" class="mean"/>')

    output = OUT_DIR / f"{gene_name}_replicate_dot_plot.svg"
    output.write_text("\n".join(parts + ["</svg>\n"]))


def create_orf_tables(annotation: pd.DataFrame, cds: pd.DataFrame, expr: pd.DataFrame) -> None:
    annotated = annotation.merge(cds, on=["transcript_id", "protein_id"], how="left")
    selected = []
    for gene_name in GENES:
        gene = prepare_gene_data(expr, annotated, gene_name)
        if gene.empty:
            continue
        keep_cols = [
            "gene_name",
            "gene_id",
            "transcript_id",
            "transcript_name",
            "transcript_type",
            "protein_id",
            "cds_bp",
            "cds_segments",
            "protein_length_aa_est",
            "max_mean_TPM",
        ]
        keep_cols += [f"mean_TPM_{minute}" for minute in TIMEPOINTS]
        out = gene[keep_cols].copy()
        out.to_csv(OUT_DIR / f"{gene_name}_orf_cds_annotation.csv", index=False)
        selected.append(out)

    all_selected = pd.concat(selected, ignore_index=True)
    length_summary = (
        all_selected[all_selected["protein_id"].fillna("").ne("")]
        .groupby("gene_name")
        .agg(
            protein_coding_transcripts=("transcript_id", "nunique"),
            distinct_estimated_protein_lengths=("protein_length_aa_est", "nunique"),
            min_protein_length_aa=("protein_length_aa_est", "min"),
            max_protein_length_aa=("protein_length_aa_est", "max"),
            max_mean_TPM=("max_mean_TPM", "max"),
        )
        .reset_index()
        .sort_values(["distinct_estimated_protein_lengths", "max_mean_TPM"], ascending=False)
    )
    length_summary.to_csv(OUT_DIR / "selected_genes_orf_length_summary.csv", index=False)


def write_readme() -> None:
    lines = [
        "# Isoform Validation Outputs",
        "",
        "## Replicate Dot Plots",
        "",
        "These plots answer:",
        "",
        "`Are the isoform TPM changes visible in individual replicates, or only in the mean?`",
        "",
        "Each dot is one replicate. The black line connects condition means.",
        "",
        "A strong candidate should usually show dots moving in the same direction across most or all replicates.",
        "",
        "## ORF/CDS Tables",
        "",
        "These tables answer:",
        "",
        "`Could the changed transcript isoform encode a different protein product?`",
        "",
        "Useful columns:",
        "",
        "- `transcript_type`: protein_coding, retained_intron, nonsense_mediated_decay, etc.",
        "- `protein_id`: Ensembl protein ID, present when a protein product is annotated.",
        "- `cds_bp`: total coding sequence length in base pairs from GENCODE CDS records.",
        "- `protein_length_aa_est`: estimated amino acid length from CDS length.",
        "",
        "If two isoforms from the same gene have different protein IDs or different estimated protein lengths, the isoform change may affect the encoded protein.",
        "",
        "Protein domain effects require an additional domain database such as UniProt, Pfam, or InterPro. The current outputs are ORF/CDS-level preliminary checks.",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expr = build_expression_table()
    annotation = load_transcript_annotation()
    cds = load_cds_annotation()

    for gene_name in GENES:
        gene = prepare_gene_data(expr, annotation, gene_name)
        if gene.empty:
            print(f"Skipping {gene_name}")
            continue
        create_replicate_dot_plot(gene_name, gene)
        print(f"Wrote replicate plot for {gene_name}")

    create_orf_tables(annotation, cds, expr)
    write_readme()
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
