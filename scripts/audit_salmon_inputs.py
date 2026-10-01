#!/usr/bin/env python3
"""Audit all deposited Salmon quantifications without biological filtering."""

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "metadata" / "samples.tsv"
OUTDIR = PROJECT / "results" / "qc"
REQUIRED_COLUMNS = {"Name", "Length", "EffectiveLength", "TPM", "NumReads"}


def main() -> None:
    samples = pd.read_csv(MANIFEST, sep="\t")
    if len(samples) != 12:
        raise ValueError(f"Expected 12 samples, found {len(samples)}")
    if samples["sample_id"].duplicated().any():
        raise ValueError("sample_id values must be unique")
    if sorted(samples.groupby("time_min").size().tolist()) != [3, 3, 3, 3]:
        raise ValueError("Expected three samples at each of 0, 12, 30 and 60 min")

    summaries = []
    reference_ids = None
    log_tpm = {}

    for row in samples.itertuples(index=False):
        path = PROJECT / row.quant_file
        if not path.exists():
            raise FileNotFoundError(path)
        quant = pd.read_csv(path, sep="\t")
        missing = REQUIRED_COLUMNS.difference(quant.columns)
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        if quant["Name"].duplicated().any():
            raise ValueError(f"{path.name}: duplicate transcript identifiers")
        if (quant[["Length", "EffectiveLength", "TPM", "NumReads"]] < 0).any().any():
            raise ValueError(f"{path.name}: negative abundance or length value")

        ids = quant["Name"].astype(str)
        if reference_ids is None:
            reference_ids = ids.tolist()
        elif ids.tolist() != reference_ids:
            raise ValueError(f"{path.name}: transcript identifiers/order differ")

        summaries.append(
            {
                "sample_id": row.sample_id,
                "geo_accession": row.geo_accession,
                "time_min": row.time_min,
                "replicate": row.replicate,
                "n_transcripts": len(quant),
                "n_nonzero_tpm": int((quant["TPM"] > 0).sum()),
                "tpm_sum": quant["TPM"].sum(),
                "estimated_assigned_fragments": quant["NumReads"].sum(),
            }
        )
        log_tpm[row.sample_id] = np.log2(quant["TPM"].to_numpy() + 0.1)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUTDIR / "salmon_input_audit.tsv", sep="\t", index=False)

    matrix = pd.DataFrame(log_tpm)
    corr = matrix.corr(method="pearson")
    corr.index.name = "sample_id"
    corr.to_csv(OUTDIR / "sample_log2tpm_correlations.tsv", sep="\t")

    status = [
        "# Salmon input audit",
        "",
        "This audit applies no expression, effect-size or significance threshold.",
        "All 12 deposited samples are retained.",
        "",
        f"- Samples: {len(summary)}",
        f"- Time points: {', '.join(map(str, sorted(summary.time_min.unique())))} min",
        f"- Replicates per time point: {summary.groupby('time_min').size().min()}",
        f"- Transcripts per sample: {summary.n_transcripts.min():,}",
        f"- TPM-sum range: {summary.tpm_sum.min():,.3f}–{summary.tpm_sum.max():,.3f}",
        "- Transcript identifiers and row order: identical across samples",
        f"- Pairwise log2(TPM + 0.1) correlation range: "
        f"{corr.where(~np.eye(len(corr), dtype=bool)).stack().min():.3f}–"
        f"{corr.where(~np.eye(len(corr), dtype=bool)).stack().max():.3f}",
        "",
        "Passing this audit does not establish biological comparability or DTU.",
        "It only establishes that the deposited quantification files are complete",
        "and structurally compatible for downstream modelling.",
    ]
    (OUTDIR / "salmon_input_audit.md").write_text("\n".join(status) + "\n")
    print("\n".join(status))


if __name__ == "__main__":
    main()
