#!/usr/bin/env python3
"""
Analyze transcript-level RNA-seq quantification from Kaya & Adebali 2025.

Question:
Which transcript isoforms change after UV exposure compared with the non-UV
control at 12, 30, and 60 minutes?

Input:
geo_rnaseq_quant/*_quant.sf.gz files downloaded from GEO accession GSE268349.

Output:
results/transcript_isoform_uv_vs_control_all.csv
results/transcript_isoform_uv_vs_control_changed.csv
results/transcript_isoform_uv_vs_control_summary.csv
"""

from pathlib import Path
import gzip
import re

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
QUANT_DIR = PROJECT_DIR / "geo_rnaseq_quant"
RESULTS_DIR = PROJECT_DIR / "results"
GENCODE_GTF = PROJECT_DIR / "gencode.v35.annotation.gtf.gz"

# Small pseudocount prevents division by zero when a transcript has TPM = 0.
PSEUDOCOUNT = 0.1

# A transcript must be expressed in at least one condition to be interpreted.
MIN_MAX_MEAN_TPM = 1.0

# abs(log2FC) >= 1 means at least a 2-fold increase or decrease.
MIN_ABS_MEAN_LOG2FC = 1.0

# We have three biological replicates. Require all three to move the same way.
MIN_CONSISTENT_REPLICATES = 3


SAMPLES = [
    ("GSM8289927_rna_0_1_quant.sf.gz", "control", 0, 1),
    ("GSM8289931_rna_0_2_quant.sf.gz", "control", 0, 2),
    ("GSM8289935_rna_0_3_quant.sf.gz", "control", 0, 3),
    ("GSM8289928_rna_12_1_quant.sf.gz", "uv_12min", 12, 1),
    ("GSM8289932_rna_12_2_quant.sf.gz", "uv_12min", 12, 2),
    ("GSM8289936_rna_12_3_quant.sf.gz", "uv_12min", 12, 3),
    ("GSM8289929_rna_30_1_quant.sf.gz", "uv_30min", 30, 1),
    ("GSM8289933_rna_30_2_quant.sf.gz", "uv_30min", 30, 2),
    ("GSM8289937_rna_30_3_quant.sf.gz", "uv_30min", 30, 3),
    ("GSM8289930_rna_60_1_quant.sf.gz", "uv_60min", 60, 1),
    ("GSM8289934_rna_60_2_quant.sf.gz", "uv_60min", 60, 2),
    ("GSM8289938_rna_60_3_quant.sf.gz", "uv_60min", 60, 3),
]

ANNOTATION_COLUMNS = [
    "transcript_id",
    "gene_id",
    "gene_name",
    "gene_type",
    "transcript_name",
    "transcript_type",
    "protein_id",
]


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    """Convert a GTF attribute string into a Python dictionary."""
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', attribute_text))


def load_transcript_annotation() -> pd.DataFrame:
    """Read GENCODE v35 transcript annotations and keep one row per isoform."""
    if not GENCODE_GTF.exists():
        print(f"Annotation file not found: {GENCODE_GTF}")
        print("Continuing without gene/protein annotation.")
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    rows = []
    with gzip.open(GENCODE_GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue

            attributes = parse_gtf_attributes(fields[8])
            rows.append(
                {
                    "transcript_id": attributes.get("transcript_id", ""),
                    "gene_id": attributes.get("gene_id", ""),
                    "gene_name": attributes.get("gene_name", ""),
                    "gene_type": attributes.get("gene_type", ""),
                    "transcript_name": attributes.get("transcript_name", ""),
                    "transcript_type": attributes.get("transcript_type", ""),
                    "protein_id": attributes.get("protein_id", ""),
                }
            )

    return pd.DataFrame(rows).drop_duplicates("transcript_id")


def read_sample(filename: str, condition: str, minute: int, replicate: int) -> pd.DataFrame:
    """Read one quant.sf.gz file and keep transcript TPM/count information."""
    path = QUANT_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")

    df = pd.read_csv(path, sep="\t")
    df = df.rename(
        columns={
            "Name": "transcript_id",
            "TPM": f"{condition}_rep{replicate}_TPM",
            "NumReads": f"{condition}_rep{replicate}_NumReads",
        }
    )
    return df[
        [
            "transcript_id",
            "Length",
            "EffectiveLength",
            f"{condition}_rep{replicate}_TPM",
            f"{condition}_rep{replicate}_NumReads",
        ]
    ]


def build_expression_table() -> pd.DataFrame:
    """Merge all 12 samples into one transcript-by-sample table."""
    merged = None
    for filename, condition, minute, replicate in SAMPLES:
        sample_df = read_sample(filename, condition, minute, replicate)
        if merged is None:
            merged = sample_df
        else:
            merged = merged.merge(
                sample_df.drop(columns=["Length", "EffectiveLength"]),
                on="transcript_id",
                how="outer",
            )

    return merged.fillna(0)


def summarize_timepoint(expr: pd.DataFrame, minute: int) -> pd.DataFrame:
    """Calculate UV-vs-control transcript changes for one time point."""
    uv_prefix = f"uv_{minute}min"
    control_tpm_cols = [f"control_rep{i}_TPM" for i in (1, 2, 3)]
    uv_tpm_cols = [f"{uv_prefix}_rep{i}_TPM" for i in (1, 2, 3)]

    out = expr[["transcript_id", "Length", "EffectiveLength"]].copy()
    out["time_min"] = minute

    out["control_mean_TPM"] = expr[control_tpm_cols].mean(axis=1)
    out["uv_mean_TPM"] = expr[uv_tpm_cols].mean(axis=1)
    out["max_mean_TPM"] = out[["control_mean_TPM", "uv_mean_TPM"]].max(axis=1)
    out["mean_log2FC"] = np.log2(
        (out["uv_mean_TPM"] + PSEUDOCOUNT)
        / (out["control_mean_TPM"] + PSEUDOCOUNT)
    )

    for replicate in (1, 2, 3):
        control_col = f"control_rep{replicate}_TPM"
        uv_col = f"{uv_prefix}_rep{replicate}_TPM"
        out[f"log2FC_rep{replicate}"] = np.log2(
            (expr[uv_col] + PSEUDOCOUNT) / (expr[control_col] + PSEUDOCOUNT)
        )

    lfc_cols = [f"log2FC_rep{i}" for i in (1, 2, 3)]
    out["sd_log2FC"] = out[lfc_cols].std(axis=1)

    out["direction"] = np.where(
        out["mean_log2FC"] > 0,
        "up",
        np.where(out["mean_log2FC"] < 0, "down", "unchanged"),
    )
    out["consistent_replicates"] = np.where(
        out["direction"] == "up",
        (out[lfc_cols] > 0).sum(axis=1),
        np.where(out["direction"] == "down", (out[lfc_cols] < 0).sum(axis=1), 0),
    )

    out["passes_expression_filter"] = out["max_mean_TPM"] >= MIN_MAX_MEAN_TPM
    out["passes_fold_change_filter"] = out["mean_log2FC"].abs() >= MIN_ABS_MEAN_LOG2FC
    out["passes_consistency_filter"] = (
        out["consistent_replicates"] >= MIN_CONSISTENT_REPLICATES
    )
    out["changed"] = (
        out["passes_expression_filter"]
        & out["passes_fold_change_filter"]
        & out["passes_consistency_filter"]
    )

    return out.sort_values("mean_log2FC", ascending=False)


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    expr = build_expression_table()
    annotation = load_transcript_annotation()
    all_results = pd.concat(
        [summarize_timepoint(expr, minute) for minute in (12, 30, 60)],
        ignore_index=True,
    )

    if not annotation.empty:
        all_results = annotation.merge(all_results, on="transcript_id", how="right")

    changed = all_results[all_results["changed"]].copy()
    protein_coding_changed = changed[
        (changed["gene_type"] == "protein_coding")
        & (changed["transcript_type"] == "protein_coding")
        & changed["protein_id"].fillna("").ne("")
    ].copy()
    summary = (
        changed.groupby(["time_min", "direction"])
        .size()
        .reset_index(name="changed_transcript_count")
        .sort_values(["time_min", "direction"])
    )
    protein_coding_summary = (
        protein_coding_changed.groupby(["time_min", "direction"])
        .size()
        .reset_index(name="changed_protein_coding_transcript_count")
        .sort_values(["time_min", "direction"])
    )

    all_path = RESULTS_DIR / "transcript_isoform_uv_vs_control_all.csv"
    changed_path = RESULTS_DIR / "transcript_isoform_uv_vs_control_changed.csv"
    summary_path = RESULTS_DIR / "transcript_isoform_uv_vs_control_summary.csv"
    protein_path = (
        RESULTS_DIR / "protein_coding_transcript_isoforms_uv_vs_control_changed.csv"
    )
    protein_summary_path = (
        RESULTS_DIR / "protein_coding_transcript_isoforms_uv_vs_control_summary.csv"
    )

    all_results.to_csv(all_path, index=False)
    changed.to_csv(changed_path, index=False)
    summary.to_csv(summary_path, index=False)
    protein_coding_changed.to_csv(protein_path, index=False)
    protein_coding_summary.to_csv(protein_summary_path, index=False)

    print(f"Wrote {all_path}")
    print(f"Wrote {changed_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {protein_path}")
    print(f"Wrote {protein_summary_path}")
    print()
    print("All changed transcript isoforms:")
    print(summary.to_string(index=False))
    print()
    print("Protein-coding changed transcript isoforms with protein IDs:")
    print(protein_coding_summary.to_string(index=False))


if __name__ == "__main__":
    main()
