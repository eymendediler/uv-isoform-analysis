#!/usr/bin/env python3
"""Compare SUPPA2 DTU calls with the existing exact-permutation candidates."""

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_DIR / "results"
SUPPA_DIR = RESULTS_DIR / "suppa2"


def main() -> None:
    suppa = pd.read_csv(SUPPA_DIR / "suppa2_isoform_dtu_candidates.csv")
    permutation = pd.read_csv(RESULTS_DIR / "statistical_DTU_candidate_transcripts.csv")

    rows = []
    overlap_frames = []
    for minute in (12, 30, 60):
        s = suppa.loc[suppa["time_min"].eq(minute)].copy()
        p = permutation.loc[permutation["time_min"].eq(minute)].copy()
        s_ids = set(s["transcript_id"])
        p_ids = set(p["transcript_id"])
        overlap = s_ids & p_ids

        joined = s.loc[s["transcript_id"].isin(overlap)].merge(
            p[["transcript_id", "delta_usage"]],
            on="transcript_id",
            how="inner",
        )
        joined["direction_agrees"] = (
            joined["delta_psi"].mul(joined["delta_usage"]).gt(0)
        )
        joined["time_min"] = minute
        overlap_frames.append(joined)

        rows.append(
            {
                "time_min": minute,
                "suppa2_candidates": len(s_ids),
                "permutation_candidates": len(p_ids),
                "overlap": len(overlap),
                "jaccard": len(overlap) / len(s_ids | p_ids),
                "direction_agreement": joined["direction_agrees"].mean(),
            }
        )

    pd.DataFrame(rows).to_csv(
        SUPPA_DIR / "method_concordance_summary.csv", index=False
    )
    pd.concat(overlap_frames, ignore_index=True).to_csv(
        SUPPA_DIR / "cross_method_supported_isoforms.csv", index=False
    )

    recurrence = (
        suppa.groupby(
            [
                "transcript_id",
                "gene_id",
                "gene_name",
                "gene_type",
                "transcript_name",
                "transcript_type",
                "protein_id",
            ],
            dropna=False,
        )
        .agg(
            timepoints=("time_min", lambda x: ",".join(map(str, sorted(set(x))))),
            n_timepoints=("time_min", "nunique"),
            max_abs_delta_psi=("delta_psi", lambda x: x.abs().max()),
            min_suppa_pvalue=("suppa_pvalue", "min"),
        )
        .reset_index()
    )
    recurrence = recurrence.loc[recurrence["n_timepoints"].ge(2)].sort_values(
        ["n_timepoints", "max_abs_delta_psi"], ascending=[False, False]
    )
    recurrence.to_csv(SUPPA_DIR / "recurrent_isoform_dtu_candidates.csv", index=False)

    print(pd.DataFrame(rows).to_string(index=False))
    print(f"Cross-method supported rows: {sum(len(x) for x in overlap_frames)}")
    print(f"Recurrent SUPPA2 isoforms: {len(recurrence)}")


if __name__ == "__main__":
    main()
