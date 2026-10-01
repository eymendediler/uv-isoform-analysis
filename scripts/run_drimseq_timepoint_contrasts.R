#!/usr/bin/env Rscript

# Secondary localization analysis for the primary categorical-time DTU model.
# The omnibus stageR analysis remains the discovery analysis. Here, each UV
# recovery time is contrasted with 0 min using the already fitted strict
# DRIMSeq model. Multiplicity correction is global across all transcript x
# time contrasts, rather than only across preselected candidate genes.

suppressPackageStartupMessages({
  library(data.table)
  library(DRIMSeq)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
project <- normalizePath(file.path(dirname(file_arg), ".."))
run_dir <- file.path(project, "results", "dtu_v2", "strict_unblocked_scaled_tpm")
fit_file <- file.path(run_dir, "drimseq_fit.rds")
usage_file <- file.path(run_dir, "transcript_omnibus_and_usage.tsv.gz")
annotation_file <- file.path(project, "metadata", "gencode_v35_transcripts.tsv.gz")
out_file <- file.path(
  project, "results", "dtu_v2",
  "timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1.tsv.gz"
)
summary_file <- file.path(
  project, "results", "dtu_v2",
  "timepoint_contrasts_strict_scaled_tpm_prop_sd_0_1_summary.tsv"
)

d <- readRDS(fit_file)
usage <- as.data.table(read.delim(gzfile(usage_file), check.names = FALSE))
annotation <- as.data.table(read.delim(gzfile(annotation_file), check.names = FALSE))

# The official rnaseqDTU workflow recommends removing transcript tests whose
# observed within-gene proportion has SD < 0.1 before multiple-testing control.
count_dt <- as.data.table(counts(d))
sample_cols <- setdiff(names(count_dt), c("gene_id", "feature_id"))
cts <- as.matrix(count_dt[, ..sample_cols])
gene_cts <- rowsum(cts, count_dt$gene_id)
totals <- gene_cts[match(count_dt$gene_id, rownames(gene_cts)), , drop = FALSE]
props <- cts / totals
prop_sd <- apply(props, 1L, sd)
prop_sd_table <- count_dt[, .(gene_id, feature_id)]
prop_sd_table[, proportion_sd_all_samples := prop_sd]

coef_map <- data.table(
  coefficient = c("time12", "time30", "time60"),
  time_min = c(12L, 30L, 60L)
)

contrast_results <- lapply(seq_len(nrow(coef_map)), function(i) {
  tested <- dmTest(d, coef = coef_map$coefficient[i], one_way = TRUE)
  tab <- as.data.table(DRIMSeq::results(tested, level = "feature"))
  tab[, `:=`(
    coefficient = coef_map$coefficient[i],
    time_min = coef_map$time_min[i]
  )]
  tab
})
contrast_results <- rbindlist(contrast_results, use.names = TRUE)
contrast_results <- merge(
  contrast_results,
  prop_sd_table,
  by = c("gene_id", "feature_id"),
  all.x = TRUE
)
contrast_results[is.na(pvalue), pvalue := 1]
contrast_results[, removed_by_prop_sd_0_1 :=
                   is.na(proportion_sd_all_samples) |
                   proportion_sd_all_samples < 0.10]
contrast_results[removed_by_prop_sd_0_1 == TRUE, pvalue := 1]

# One correction family: every retained transcript at all three time points.
# This is deliberately more conservative than correcting selected genes only.
contrast_results[, global_bh_qvalue := p.adjust(pvalue, method = "BH")]

usage_cols <- c(
  "gene_id", "feature_id", "tpm_mean_usage_0",
  paste0("tpm_mean_usage_", c(12L, 30L, 60L)),
  paste0("tpm_delta_usage_", c(12L, 30L, 60L), "_vs_0")
)
contrast_results <- merge(
  contrast_results,
  usage[, ..usage_cols],
  by = c("gene_id", "feature_id"),
  all.x = TRUE
)
contrast_results[, tpm_mean_usage_time := fifelse(
  time_min == 12L, tpm_mean_usage_12,
  fifelse(time_min == 30L, tpm_mean_usage_30, tpm_mean_usage_60)
)]
contrast_results[, tpm_delta_usage_vs_0 := fifelse(
  time_min == 12L, tpm_delta_usage_12_vs_0,
  fifelse(time_min == 30L, tpm_delta_usage_30_vs_0,
          tpm_delta_usage_60_vs_0)
)]
contrast_results[, direction := fifelse(
  tpm_delta_usage_vs_0 > 0, "increased usage",
  fifelse(tpm_delta_usage_vs_0 < 0, "decreased usage", "no change")
)]

contrast_results <- merge(
  contrast_results,
  annotation[, .(
    feature_id = transcript_id, gene_name, transcript_name,
    gene_type, transcript_type, protein_id
  )],
  by = "feature_id",
  all.x = TRUE
)
setorder(contrast_results, global_bh_qvalue, pvalue)
fwrite(contrast_results, out_file, sep = "\t")

summary <- contrast_results[, .(
  tested_transcript_time_contrasts = .N,
  unique_transcripts = uniqueN(feature_id),
  unique_genes = uniqueN(gene_id),
  removed_transcript_time_contrasts_by_prop_sd = sum(removed_by_prop_sd_0_1),
  global_bh_q_le_0_05 = sum(global_bh_qvalue <= 0.05),
  genes_with_global_bh_q_le_0_05 = uniqueN(gene_id[global_bh_qvalue <= 0.05])
)]
fwrite(summary, summary_file, sep = "\t")
print(summary)
