#!/usr/bin/env python3
"""Run SUPPA2 transcript-usage and local-splicing analyses on the Salmon data.

SUPPA2 itself is intentionally not vendored. Clone the official repository to
``.tools/SUPPA`` (or set ``SUPPA_DIR``), then run this script from the project
environment. Intermediate files remain under ``results/suppa2_work`` and are
ignored by Git; compact candidate and summary tables are committed.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_transcript_isoforms import SAMPLES, load_transcript_annotation


PROJECT_DIR = Path(__file__).resolve().parents[1]
QUANT_DIR = PROJECT_DIR / "geo_rnaseq_quant"
GTF_GZ = PROJECT_DIR / "gencode.v35.annotation.gtf.gz"
RESULTS_DIR = PROJECT_DIR / "results" / "suppa2"
WORK_DIR = PROJECT_DIR / "results" / "suppa2_work"
SUPPA_DIR = Path(os.environ.get("SUPPA_DIR", PROJECT_DIR / ".tools" / "SUPPA"))
SUPPA = SUPPA_DIR / "suppa.py"

TIMEPOINTS = (12, 30, 60)
EVENT_TYPES = ("SE", "A3", "A5", "MX", "RI", "AF", "AL")
MIN_ABS_DPSI = 0.10
MAX_PVALUE = 0.05


def run(*args: str) -> None:
    """Run SUPPA2 and fail with the exact command if it exits unsuccessfully."""
    command = [sys.executable, str(SUPPA), *map(str, args)]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def apply_performance_patch() -> None:
    """Apply a semantics-preserving fix for SUPPA2's quadratic nearest lookup."""
    diff_tools = SUPPA_DIR / "lib" / "diff_tools.py"
    source = diff_tools.read_text()
    replacement = (
        "bisect_left(replicates_logtpms, close_rep_logtpm)"
    )
    if replacement in source:
        return
    original = "replicates_logtpms.index(close_rep_logtpm)"
    if original not in source:
        raise RuntimeError(
            "SUPPA2 lookup code differs from the supported v2.4 source; "
            "refusing to patch an unknown implementation."
        )
    diff_tools.write_text(
        source.replace(original, replacement, 1)
    )


def write_expression_files() -> dict[str, Path]:
    """Create one SUPPA-compatible TPM matrix per condition."""
    per_condition: dict[str, list[tuple[str, pd.Series]]] = {}
    for filename, condition, _minute, replicate in SAMPLES:
        sample_name = f"{condition}_rep{replicate}"
        quant = pd.read_csv(QUANT_DIR / filename, sep="\t", compression="gzip")
        per_condition.setdefault(condition, []).append(
            (sample_name, quant.set_index("Name")["TPM"])
        )

    paths = {}
    for condition, samples in per_condition.items():
        matrix = pd.concat([series.rename(name) for name, series in samples], axis=1)
        path = WORK_DIR / f"{condition}.tpm"
        matrix.to_csv(path, sep="\t", index=True, index_label=False)
        paths[condition] = path
    return paths


def split_psi(all_psi: Path, expression_paths: dict[str, Path], stem: str) -> dict[str, Path]:
    """Split a 12-sample SUPPA PSI matrix into condition-specific matrices."""
    psi = pd.read_csv(all_psi, sep="\t", index_col=0)
    paths = {}
    for condition, expression_path in expression_paths.items():
        sample_names = list(pd.read_csv(expression_path, sep="\t", nrows=0).columns)
        path = WORK_DIR / f"{stem}_{condition}.psi"
        psi[sample_names].to_csv(path, sep="\t", index=True, index_label=False)
        paths[condition] = path
    return paths


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg correction that preserves missing values."""
    output = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    p = values.loc[valid].astype(float).to_numpy()
    if not len(p):
        return output
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    restored = np.empty(len(p))
    restored[order] = np.clip(adjusted, 0, 1)
    output.loc[valid] = restored
    return output


def collect_dpsi(prefix: Path, minute: int, analysis: str) -> pd.DataFrame:
    """Read a SUPPA dPSI file and standardize its columns."""
    path = prefix.with_suffix(".dpsi")
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.columns = ["delta_psi", "suppa_pvalue"]
    frame.index.name = "event_id"
    frame = frame.reset_index()
    frame["time_min"] = minute
    frame["analysis"] = analysis
    frame["padj_BH"] = bh_adjust(frame["suppa_pvalue"])
    frame["candidate"] = (
        frame["delta_psi"].abs().ge(MIN_ABS_DPSI)
        & frame["suppa_pvalue"].le(MAX_PVALUE)
    )
    return frame


def annotate_isoforms(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach GENCODE metadata to SUPPA transcript events."""
    frame["transcript_id"] = frame["event_id"].str.split(";", n=1).str[-1]
    annotation = load_transcript_annotation()
    return annotation.merge(frame, on="transcript_id", how="right")


def summarize(candidates: pd.DataFrame, label: str) -> pd.DataFrame:
    """Count candidates by time point and event type."""
    group_cols = ["time_min"]
    if label == "local_event":
        candidates["event_type"] = candidates["event_id"].str.extract(r";([A-Z0-9]+):")
        group_cols.append("event_type")
    summary = (
        candidates.groupby(group_cols, dropna=False)
        .size()
        .reset_index(name="candidate_count")
        .sort_values(group_cols)
    )
    summary.insert(0, "analysis", label)
    return summary


def run_comparisons(
    iox: Path,
    psi_paths: dict[str, Path],
    expression_paths: dict[str, Path],
    analysis: str,
) -> pd.DataFrame:
    """Compare each UV time point with the shared control."""
    outputs = []
    for minute in TIMEPOINTS:
        uv = f"uv_{minute}min"
        prefix = WORK_DIR / f"{analysis}_control_vs_{minute}min"
        run(
            "diffSplice",
            "--method", "empirical",
            "--input", iox,
            "--psi", psi_paths["control"], psi_paths[uv],
            "--tpm", expression_paths["control"], expression_paths[uv],
            "--area", "1000",
            "--lower-bound", "0.05",
            "--tpm-threshold", "1",
            "--output", prefix,
        )
        outputs.append(collect_dpsi(prefix, minute, analysis))
    return pd.concat(outputs, ignore_index=True)


def main() -> None:
    if not SUPPA.exists():
        raise FileNotFoundError(
            f"SUPPA2 not found at {SUPPA}. Clone https://github.com/comprna/SUPPA "
            "there or set SUPPA_DIR."
        )
    if not GTF_GZ.exists():
        raise FileNotFoundError(f"Missing annotation: {GTF_GZ}")
    apply_performance_patch()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    gtf = WORK_DIR / "gencode.v35.annotation.gtf"
    with gzip.open(GTF_GZ, "rb") as source, gtf.open("wb") as destination:
        shutil.copyfileobj(source, destination)

    expression_paths = write_expression_files()
    all_expression = WORK_DIR / "all_samples.tpm"
    matrices = [
        pd.read_csv(path, sep="\t", index_col=0) for path in expression_paths.values()
    ]
    pd.concat(matrices, axis=1).to_csv(
        all_expression, sep="\t", index=True, index_label=False
    )

    iso_ioi_prefix = WORK_DIR / "gencode_isoforms"
    run("generateEvents", "--input-file", gtf, "--output-file", iso_ioi_prefix, "--format", "ioi")
    iso_ioi = iso_ioi_prefix.with_suffix(".ioi")
    iso_psi_prefix = WORK_DIR / "isoform_all"
    run("psiPerIsoform", "--gtf-file", gtf, "--expression-file", all_expression, "--output-file", iso_psi_prefix)
    iso_psi_paths = split_psi(
        Path(f"{iso_psi_prefix}_isoform.psi"), expression_paths, "isoform"
    )
    isoforms = annotate_isoforms(
        run_comparisons(iso_ioi, iso_psi_paths, expression_paths, "isoform")
    )
    isoform_candidates = isoforms.loc[isoforms["candidate"]].copy()
    isoform_candidates.to_csv(RESULTS_DIR / "suppa2_isoform_dtu_candidates.csv", index=False)

    event_prefix = WORK_DIR / "gencode_events"
    run(
        "generateEvents", "--input-file", gtf, "--output-file", event_prefix,
        "--format", "ioe", "--event-type", "SE", "SS", "MX", "RI", "FL",
    )
    event_frames = []
    for event_type in EVENT_TYPES:
        ioe = Path(f"{event_prefix}_{event_type}_strict.ioe")
        if not ioe.exists():
            continue
        psi_prefix = WORK_DIR / f"event_{event_type}_all"
        run("psiPerEvent", "--ioe-file", ioe, "--expression-file", all_expression, "--output-file", psi_prefix)
        psi_paths = split_psi(psi_prefix.with_suffix(".psi"), expression_paths, f"event_{event_type}")
        event_frames.append(run_comparisons(ioe, psi_paths, expression_paths, f"local_{event_type}"))

    events = pd.concat(event_frames, ignore_index=True)
    event_candidates = events.loc[events["candidate"]].copy()
    event_candidates.to_csv(RESULTS_DIR / "suppa2_local_splicing_candidates.csv", index=False)

    summary = pd.concat(
        [
            summarize(isoform_candidates, "isoform"),
            summarize(event_candidates, "local_event"),
        ],
        ignore_index=True,
    )
    summary.to_csv(RESULTS_DIR / "suppa2_candidate_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
