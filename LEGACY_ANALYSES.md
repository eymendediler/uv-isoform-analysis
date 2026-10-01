# Superseded preliminary analyses

The files listed here are retained for provenance, auditability and historical
reproduction. They are **not the current primary results** and should not be
used for biological conclusions or presentations without an explicit legacy
label.

## Why they were superseded

The earliest analysis used pairwise label permutations on descriptive usage
quantities and later compared them with SUPPA2. It was useful for initial
candidate generation, but it did not provide the preferred count-based
replicate model for the four-time-point design. Some earlier candidate lists
also preceded the recommended condition-blind DRIMSeq proportion-SD filter.

The current workflow instead uses all replicate-level samples in a categorical
DRIMSeq Dirichlet-multinomial time model, stageR hierarchical error control,
filter-sensitivity analysis, protein-sequence validation, and globally adjusted
time-point contrasts.

## Legacy locations

- Top-level `results/statistical_DTE_*` and `results/statistical_DTU_*` tables
- Top-level `results/*isoform_uv_vs_control*` tables
- `results/suppa2/`
- `preliminary_data/`
- `reproducibility_checks/` relating to the preliminary/SUPPA2 workflow
- Scripts identified as `legacy` in `scripts/README.md`

These files have not been deleted because preserving unsuccessful and
superseded approaches prevents accidental rewriting of the analysis history.
The authoritative entry point is `results/dtu_v2/analysis_status.md`.
