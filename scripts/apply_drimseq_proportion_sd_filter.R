#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DRIMSeq)
  library(stageR)
})

project <- normalizePath(getwd(), mustWork = TRUE)
root <- file.path(project, "results", "dtu_v2")
runs <- paste0(c("loose", "medium", "strict"), "_unblocked_scaled_tpm")
alpha <- 0.05
sd_cutoff <- 0.10

for (run in runs) {
  folder <- file.path(root, run)
  d <- readRDS(file.path(folder, "drimseq_fit.rds"))
  gene_res <- as.data.table(DRIMSeq::results(d, level = "gene"))
  tx_res <- as.data.table(DRIMSeq::results(d, level = "feature"))
  gene_res[is.na(pvalue), pvalue := 1]
  tx_res[is.na(pvalue), pvalue := 1]

  count_dt <- as.data.table(counts(d))
  sample_cols <- setdiff(names(count_dt), c("gene_id", "feature_id"))
  cts <- as.matrix(count_dt[, ..sample_cols])
  gene_cts <- rowsum(cts, count_dt$gene_id)
  totals <- gene_cts[match(count_dt$gene_id, rownames(gene_cts)), , drop = FALSE]
  props <- cts / totals
  prop_sd <- apply(props, 1L, sd)
  low_sd <- prop_sd < sd_cutoff | is.na(prop_sd)

  tx_res[, proportion_sd_all_samples := prop_sd]
  tx_res[, removed_by_prop_sd_0_1 := low_sd]
  tx_res[low_sd, pvalue := 1]

  p_screen <- gene_res$pvalue
  names(p_screen) <- gene_res$gene_id
  p_confirmation <- matrix(
    tx_res$pvalue,
    ncol = 1L,
    dimnames = list(tx_res$feature_id, "omnibus_time")
  )
  tx2gene <- as.data.frame(tx_res[, .(feature_id, gene_id)])
  stage <- stageRTx(
    pScreen = p_screen,
    pConfirmation = p_confirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  )
  stage <- stageWiseAdjustment(stage, method = "dtu", alpha = alpha)
  adjusted <- as.data.table(getAdjustedPValues(
    stage, onlySignificantGenes = FALSE, order = FALSE
  ))
  adjusted <- merge(
    adjusted,
    tx_res[, .(txID = feature_id, proportion_sd_all_samples, removed_by_prop_sd_0_1)],
    by = "txID",
    all.x = TRUE
  )
  fwrite(
    adjusted,
    file.path(folder, "stageR_adjusted_prop_sd_0_1.tsv"),
    sep = "\t"
  )
  summary <- data.table(
    run = run,
    tested_transcripts = nrow(tx_res),
    removed_by_prop_sd_0_1 = sum(low_sd),
    stageR_confirmed_transcripts_after_filter = sum(
      adjusted$transcript <= alpha, na.rm = TRUE
    ),
    stageR_confirmed_genes_after_filter = uniqueN(
      adjusted[transcript <= alpha, geneID]
    )
  )
  fwrite(summary, file.path(folder, "prop_sd_filter_summary.tsv"), sep = "\t")
  print(summary)
}
