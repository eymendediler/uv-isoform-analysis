#!/usr/bin/env python3
"""
Statistical transcript isoform analysis for UV time-course RNA-seq data.

This script performs two complementary analyses:

1. DTE: Differential transcript expression
   Does an isoform's expression level change after UV?

2. DTU: Differential transcript usage
   Within the same gene, does the isoform's proportion change after UV?

Because this teaching environment has only pandas/numpy available and each
condition has 3 replicates, we use an exact label-permutation test. This is
transparent and easy to understand, but a publication-level analysis should
also be repeated with specialist RNA-seq tools such as DRIMSeq/DEXSeq/satuRn
or tximport + DESeq2/edgeR where possible.
"""

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_transcript_isoforms import (
    PSEUDOCOUNT,
    build_expression_table,
    load_transcript_annotation,
)


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"

MIN_MAX_MEAN_TPM = 2.0
MIN_ABS_LOG2FC = 1.0
MIN_ABS_DELTA_USAGE = 0.10
MIN_GENE_MEAN_TPM = 2.0
MAX_PERMUTATION_P = 0.10

CONTROL_COLS = [f"control_rep{i}_TPM" for i in (1, 2, 3)]
TIMEPOINTS = [12, 30, 60]


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Calculate BH-adjusted p-values from a vector of p-values."""
    p = pvalues.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=pvalues.index)


def exact_permutation_pvalue(values_a: np.ndarray, values_b: np.ndarray) -> np.ndarray:
    """
    Exact two-sided permutation p-value for 3 control vs 3 UV replicates.

    For each transcript, we pool the 6 values, try every possible way to split
    them into two groups of 3, and ask how many splits create a difference at
    least as strong as the real control-vs-UV difference.
    """
    observed = np.abs(values_b.mean(axis=1) - values_a.mean(axis=1))
    pooled = np.concatenate([values_a, values_b], axis=1)
    possible_control_indices = list(combinations(range(6), 3))

    extreme_counts = np.zeros(len(observed), dtype=int)
    for control_idx in possible_control_indices:
        control_idx = np.array(control_idx)
        uv_idx = np.array([i for i in range(6) if i not in set(control_idx)])
        permuted_diff = np.abs(
            pooled[:, uv_idx].mean(axis=1) - pooled[:, control_idx].mean(axis=1)
        )
        extreme_counts += permuted_diff >= observed

    return extreme_counts / len(possible_control_indices)


def add_annotation(results: pd.DataFrame) -> pd.DataFrame:
    """Attach gene, transcript, and protein IDs from GENCODE v35."""
    annotation = load_transcript_annotation()
    if annotation.empty:
        return results
    return annotation.merge(results, on="transcript_id", how="right")


def run_differential_transcript_expression(expr: pd.DataFrame, minute: int) -> pd.DataFrame:
    """Analyze whether transcript isoform expression changes after UV."""
    uv_cols = [f"uv_{minute}min_rep{i}_TPM" for i in (1, 2, 3)]
    control = expr[CONTROL_COLS].to_numpy(dtype=float)
    uv = expr[uv_cols].to_numpy(dtype=float)

    control_log = np.log2(control + PSEUDOCOUNT)
    uv_log = np.log2(uv + PSEUDOCOUNT)

    out = expr[["transcript_id", "Length", "EffectiveLength"]].copy()
    out["time_min"] = minute
    out["control_mean_TPM"] = control.mean(axis=1)
    out["uv_mean_TPM"] = uv.mean(axis=1)
    out["max_mean_TPM"] = out[["control_mean_TPM", "uv_mean_TPM"]].max(axis=1)
    out["log2FC"] = np.log2(
        (out["uv_mean_TPM"] + PSEUDOCOUNT)
        / (out["control_mean_TPM"] + PSEUDOCOUNT)
    )
    out["direction"] = np.where(out["log2FC"] > 0, "up", "down")
    out["permutation_pvalue"] = exact_permutation_pvalue(control_log, uv_log)

    for replicate in (1, 2, 3):
        out[f"control_rep{replicate}_TPM"] = expr[f"control_rep{replicate}_TPM"]
        out[f"uv_rep{replicate}_TPM"] = expr[f"uv_{minute}min_rep{replicate}_TPM"]

    tested = out["max_mean_TPM"] >= MIN_MAX_MEAN_TPM
    out["padj_BH"] = np.nan
    out.loc[tested, "padj_BH"] = benjamini_hochberg(
        out.loc[tested, "permutation_pvalue"]
    )
    out["candidate_DTE"] = (
        tested
        & (out["log2FC"].abs() >= MIN_ABS_LOG2FC)
        & (out["permutation_pvalue"] <= MAX_PERMUTATION_P)
    )

    return out


def add_isoform_usage(expr: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Calculate isoform usage fraction inside each gene for each sample."""
    annotated = annotation[["transcript_id", "gene_id"]].merge(
        expr, on="transcript_id", how="right"
    )
    tpm_cols = [c for c in annotated.columns if c.endswith("_TPM")]
    gene_totals = annotated.groupby("gene_id", dropna=False)[tpm_cols].transform("sum")
    usage = annotated[["transcript_id", "gene_id", "Length", "EffectiveLength"]].copy()
    usage_cols = [c.replace("_TPM", "_usage") for c in tpm_cols]
    usage[usage_cols] = annotated[tpm_cols].div(gene_totals.replace(0, np.nan)).fillna(0)
    usage[[c.replace("_TPM", "_gene_TPM") for c in tpm_cols]] = gene_totals
    usage[tpm_cols] = annotated[tpm_cols]
    return usage


def run_differential_transcript_usage(
    usage: pd.DataFrame, annotation: pd.DataFrame, minute: int
) -> pd.DataFrame:
    """Analyze whether an isoform's proportion within its gene changes after UV."""
    control_usage_cols = [c.replace("_TPM", "_usage") for c in CONTROL_COLS]
    uv_tpm_cols = [f"uv_{minute}min_rep{i}_TPM" for i in (1, 2, 3)]
    uv_usage_cols = [c.replace("_TPM", "_usage") for c in uv_tpm_cols]
    control_gene_cols = [c.replace("_TPM", "_gene_TPM") for c in CONTROL_COLS]
    uv_gene_cols = [c.replace("_TPM", "_gene_TPM") for c in uv_tpm_cols]

    control_usage = usage[control_usage_cols].to_numpy(dtype=float)
    uv_usage = usage[uv_usage_cols].to_numpy(dtype=float)

    out = usage[["transcript_id", "gene_id", "Length", "EffectiveLength"]].copy()
    out["time_min"] = minute
    out["control_mean_usage"] = control_usage.mean(axis=1)
    out["uv_mean_usage"] = uv_usage.mean(axis=1)
    out["delta_usage"] = out["uv_mean_usage"] - out["control_mean_usage"]
    out["direction"] = np.where(out["delta_usage"] > 0, "usage_up", "usage_down")
    out["control_gene_mean_TPM"] = usage[control_gene_cols].mean(axis=1)
    out["uv_gene_mean_TPM"] = usage[uv_gene_cols].mean(axis=1)
    out["max_gene_mean_TPM"] = out[["control_gene_mean_TPM", "uv_gene_mean_TPM"]].max(
        axis=1
    )
    out["permutation_pvalue"] = exact_permutation_pvalue(control_usage, uv_usage)

    gene_isoform_counts = annotation.groupby("gene_id")["transcript_id"].nunique()
    out["annotated_isoforms_per_gene"] = out["gene_id"].map(gene_isoform_counts)

    tested = (
        (out["max_gene_mean_TPM"] >= MIN_GENE_MEAN_TPM)
        & (out["annotated_isoforms_per_gene"] > 1)
    )
    out["padj_BH"] = np.nan
    out.loc[tested, "padj_BH"] = benjamini_hochberg(
        out.loc[tested, "permutation_pvalue"]
    )
    out["candidate_DTU"] = (
        tested
        & (out["delta_usage"].abs() >= MIN_ABS_DELTA_USAGE)
        & (out["permutation_pvalue"] <= MAX_PERMUTATION_P)
    )

    return out


def protein_coding_subset(results: pd.DataFrame) -> pd.DataFrame:
    """Keep protein-coding transcript isoforms with protein IDs."""
    return results[
        (results["gene_type"] == "protein_coding")
        & (results["transcript_type"] == "protein_coding")
        & results["protein_id"].fillna("").ne("")
    ].copy()


def write_summary(df: pd.DataFrame, flag_col: str, direction_col: str, name: str) -> None:
    """Print and save count summary for candidate rows."""
    summary = (
        df[df[flag_col]]
        .groupby(["time_min", direction_col])
        .size()
        .reset_index(name="count")
        .sort_values(["time_min", direction_col])
    )
    path = RESULTS_DIR / f"{name}_summary.csv"
    summary.to_csv(path, index=False)
    print(f"\n{name}")
    print(summary.to_string(index=False))


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    expr = build_expression_table()
    annotation = load_transcript_annotation()

    dte = pd.concat(
        [run_differential_transcript_expression(expr, minute) for minute in TIMEPOINTS],
        ignore_index=True,
    )
    dte = add_annotation(dte)
    dte_candidates = dte[dte["candidate_DTE"]].copy()
    dte_protein_candidates = protein_coding_subset(dte_candidates)

    usage = add_isoform_usage(expr, annotation)
    dtu = pd.concat(
        [run_differential_transcript_usage(usage, annotation, minute) for minute in TIMEPOINTS],
        ignore_index=True,
    )
    dtu = add_annotation(dtu.drop(columns=["gene_id"]))
    dtu_candidates = dtu[dtu["candidate_DTU"]].copy()
    dtu_protein_candidates = protein_coding_subset(dtu_candidates)

    paths = {
        "dte_all": RESULTS_DIR / "statistical_DTE_all_transcripts.csv",
        "dte_candidates": RESULTS_DIR / "statistical_DTE_candidate_transcripts.csv",
        "dte_protein": RESULTS_DIR / "statistical_DTE_candidate_protein_coding_transcripts.csv",
        "dtu_all": RESULTS_DIR / "statistical_DTU_all_transcripts.csv",
        "dtu_candidates": RESULTS_DIR / "statistical_DTU_candidate_transcripts.csv",
        "dtu_protein": RESULTS_DIR / "statistical_DTU_candidate_protein_coding_transcripts.csv",
    }

    dte.to_csv(paths["dte_all"], index=False)
    dte_candidates.to_csv(paths["dte_candidates"], index=False)
    dte_protein_candidates.to_csv(paths["dte_protein"], index=False)
    dtu.to_csv(paths["dtu_all"], index=False)
    dtu_candidates.to_csv(paths["dtu_candidates"], index=False)
    dtu_protein_candidates.to_csv(paths["dtu_protein"], index=False)

    for label, path in paths.items():
        print(f"Wrote {label}: {path}")

    write_summary(dte_candidates, "candidate_DTE", "direction", "statistical_DTE_candidates")
    write_summary(
        dte_protein_candidates,
        "candidate_DTE",
        "direction",
        "statistical_DTE_protein_coding_candidates",
    )
    write_summary(dtu_candidates, "candidate_DTU", "direction", "statistical_DTU_candidates")
    write_summary(
        dtu_protein_candidates,
        "candidate_DTU",
        "direction",
        "statistical_DTU_protein_coding_candidates",
    )

    print("\nNote:")
    print("With 3 vs 3 samples, an exact permutation test has only 20 label splits.")
    print("The smallest two-sided p-value is therefore limited, so BH FDR is conservative.")
    print("Use these as statistically filtered candidates, then validate with specialist tools.")


if __name__ == "__main__":
    main()
