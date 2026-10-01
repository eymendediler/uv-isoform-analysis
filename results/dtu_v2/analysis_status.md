# DTU v2 analysis status

> **Authoritative current analysis.** Earlier top-level result CSV files,
> `results/suppa2/`, and `preliminary_data/` are superseded and retained only
> for provenance. See `../../LEGACY_ANALYSES.md`.

## Time-point localization added

The primary discovery analysis remains the categorical-time omnibus DRIMSeq +
stageR analysis. A secondary localization analysis now tests `12 vs 0`, `30 vs
0`, and `60 vs 0` using the saved strict-filter DRIMSeq fit. Transcript tests
with all-sample proportion SD below 0.10 are set to p=1, and BH correction is
applied once across all 75,477 transcript-by-time tests. This produced 316
time-localized comparisons at global BH q <= 0.05 across 158 genes.

Outputs:

- `timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1.tsv.gz`
- `primary_candidate_timepoint_statistics/primary_candidate_timepoint_heatmap.png`
- `primary_candidate_timepoint_statistics/primary_candidate_timepoint_statistics.tsv`
- replicate-aware per-gene figures in `primary_candidate_timepoint_statistics/`

These contrasts are secondary analyses and do not replace omnibus stageR
confirmation. Plot means are descriptive only; all replicate-level samples are
retained in the statistical model.

## Completed and checked

- All 12 deposited Salmon quantifications were retained.
- GENCODE v35 transcript identifiers match all 229,580 Salmon rows.
- The primary DRIMSeq input is now tximport `scaledTPM`, which preserves the
  Salmon library-size count scale while removing transcript-length bias from
  the modeled within-gene proportions. Rounded Salmon `NumReads` results are
  retained only as a sensitivity comparison and are not the primary result.
- Loose, medium and strict information filters were run independently.
- Filter retention was reproduced by an independent implementation using the
  exact DRIMSeq operation order.
- Unblocked categorical-time omnibus models completed for all three filters.
- A complementary medium-filter linear-time model used the actual ordered
  minute values 0, 12, 30 and 60. It tests a single monotonic time trend and
  does not replace the categorical omnibus model.
- Following the official Bioconductor `rnaseqDTU` workflow, transcript tests
  with across-sample proportion SD below 0.10 are now set to p=1 before stageR.
  stageR adjustment was then rerun at target OFDR 0.05 for each complete model.
- Usage effect sizes distinguish DRIMSeq model-count composition from TPM
  usage calculated against all annotated transcripts in each gene.
- In the corrected primary scaledTPM analysis, 187 transcripts were
  stageR-confirmed in all three information filters; 111 are protein-coding.
  GENCODE v35 protein sequences were compared by SHA-256: 12 genes have at
  least two robust coding transcripts that encode genuinely different
  amino-acid sequences. The earlier unfiltered values (470, 245 and 31) are
  retained only as sensitivity outputs and are not the primary result.

## Not completed

The `~ replicate + time` sensitivity model failed during DRIMSeq optimization
with a non-finite Hessian and produced no accepted result. It must not be
reported as evidence for or against DTU. A blocked analysis should be retried
with another suitable method after the experimental meaning of replicate IDs
has been confirmed with the original laboratory.

The linear-time model has so far been run only with the medium filter. Its
precision shrinkage estimate was zero, and its candidate set was substantially
larger than the categorical model's set. It is therefore retained as a
secondary exploratory analysis pending additional calibration and method
comparison, not promoted to the primary result.

## Interpretation boundary

Distinct GENCODE protein IDs indicate the potential for different encoded
protein products. Short-read RNA-seq does not prove that those protein isoforms
are translated or change at the protein level.

For example, RUNX1 has two robust transcripts encoding different protein
sequences and a maximum observed TPM-usage change of 0.521. NR4A1 has three
coding transcripts in the unfiltered sensitivity analysis, but none pass the
recommended proportion-SD filter in all three filter grids. NR4A1 should
therefore not be called a primary protein-isoform DTU candidate. RUNX1 remains
an RNA-level candidate: translation of its predicted protein isoforms is not
established.

No biological effect-size cutoff has been declared as a discovery threshold.
The tables report continuous maximum absolute TPM-usage change, with 0.05,
0.10 and 0.20 columns supplied only for transparent sensitivity summaries.
