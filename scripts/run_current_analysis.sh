#!/usr/bin/env bash

# Reproduce the accepted categorical-time DTU workflow from the repository root.
# Usage: bash scripts/run_current_analysis.sh [conda-environment-prefix]

set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
env_prefix="${1:-${project_dir}/.envs/dtu}"

cd "${project_dir}"

python scripts/audit_salmon_inputs.py
python scripts/export_gencode_annotation.py

for information_filter in loose medium strict; do
  conda run --prefix "${env_prefix}" Rscript scripts/run_drimseq_omnibus.R \
    --filter "${information_filter}" \
    --model unblocked \
    --input-scale scaled_tpm \
    --alpha 0.05 \
    --threads 1
done

conda run --prefix "${env_prefix}" Rscript \
  scripts/apply_drimseq_proportion_sd_filter.R

python scripts/summarize_dtu_sensitivity.py \
  --input-scale scaled_tpm \
  --prop-sd-filter

python scripts/validate_protein_sequences.py \
  --input-scale scaled_tpm \
  --prop-sd-filter

conda run --prefix "${env_prefix}" Rscript \
  scripts/run_gene_expression_timecourse.R

python scripts/summarize_biological_evidence.py

conda run --prefix "${env_prefix}" Rscript \
  scripts/run_drimseq_timepoint_contrasts.R

python scripts/plot_primary_dtu_statistics.py

echo "Current analysis completed successfully."
