#!/usr/bin/env Rscript

# Omnibus time-course DTU analysis. Runs one explicitly named filter/model
# combination so sensitivity analyses remain separate and auditable.

suppressPackageStartupMessages({
  library(DRIMSeq)
  library(stageR)
  library(tximport)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--filter", type = "character", default = "loose"),
  make_option("--model", type = "character", default = "unblocked"),
  make_option("--input-scale", dest = "input_scale", type = "character", default = "scaled_tpm"),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--threads", type = "integer", default = 1L)
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(opt$filter %in% c("loose", "medium", "strict"))
stopifnot(opt$model %in% c("unblocked", "blocked", "linear_time", "blocked_linear_time"))
stopifnot(opt$input_scale %in% c("scaled_tpm", "raw_counts"))
stopifnot(opt$alpha > 0, opt$alpha < 1, opt$threads >= 1L)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
project <- normalizePath(file.path(dirname(file_arg), ".."))
sample_file <- file.path(project, "metadata", "samples.tsv")
annotation_file <- file.path(project, "metadata", "gencode_v35_transcripts.tsv.gz")
run_name <- if (opt$input_scale == "raw_counts") {
  paste(opt$filter, opt$model, sep = "_")
} else {
  paste(opt$filter, opt$model, "scaled_tpm", sep = "_")
}
out_dir <- file.path(project, "results", "dtu_v2", run_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

filter_grid <- data.table(
  filter = c("loose", "medium", "strict"),
  min_samps_gene_expr = c(3L, 3L, 3L),
  min_gene_expr = c(10L, 20L, 50L),
  min_samps_feature_expr = c(2L, 3L, 3L),
  min_feature_expr = c(5L, 10L, 20L),
  min_samps_feature_prop = c(2L, 3L, 3L),
  min_feature_prop = c(0.01, 0.05, 0.10)
)
f <- filter_grid[filter == opt$filter]

samples <- fread(sample_file)
stopifnot(nrow(samples) == 12L, all(table(samples$time_min) == 3L))
annotation <- as.data.table(read.delim(gzfile(annotation_file), check.names = FALSE))

if (opt$input_scale == "scaled_tpm") {
  quant_files <- file.path(project, samples$quant_file)
  names(quant_files) <- samples$sample_id
  txi <- tximport(
    quant_files,
    type = "salmon",
    txOut = TRUE,
    countsFromAbundance = "scaledTPM",
    dropInfReps = TRUE
  )
  count_table <- as.data.table(txi$counts, keep.rownames = "feature_id")
} else {
  count_list <- lapply(seq_len(nrow(samples)), function(i) {
    q <- as.data.table(read.delim(
      gzfile(file.path(project, samples$quant_file[i])),
      check.names = FALSE
    ))[, .(Name, NumReads)]
    setnames(q, c("feature_id", samples$sample_id[i]))
    q
  })
  count_table <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), count_list)
  count_table[is.na(count_table)] <- 0
  for (column in samples$sample_id) {
    set(count_table, j = column, value = as.integer(round(count_table[[column]])))
  }
}
count_table <- merge(
  annotation[, .(feature_id = transcript_id, gene_id)],
  count_table,
  by = "feature_id",
  all.y = TRUE
)
stopifnot(!anyNA(count_table$gene_id), nrow(count_table) == 229580L)
setcolorder(count_table, c("gene_id", "feature_id", samples$sample_id))

sample_data <- data.frame(
  sample_id = samples$sample_id,
  time = factor(samples$time_min, levels = c(0, 12, 30, 60)),
  time_numeric = samples$time_min / 60,
  replicate = factor(samples$replicate),
  row.names = samples$sample_id,
  check.names = FALSE
)
d <- dmDSdata(counts = as.data.frame(count_table), samples = sample_data)
d <- dmFilter(
  d,
  min_samps_gene_expr = f$min_samps_gene_expr,
  min_gene_expr = f$min_gene_expr,
  min_samps_feature_expr = f$min_samps_feature_expr,
  min_feature_expr = f$min_feature_expr,
  min_samps_feature_prop = f$min_samps_feature_prop,
  min_feature_prop = f$min_feature_prop
)

if (opt$model == "unblocked") {
  design_full <- model.matrix(~ time, data = sample_data)
  design_null <- model.matrix(~ 1, data = sample_data)
  one_way <- TRUE
} else if (opt$model == "blocked") {
  design_full <- model.matrix(~ replicate + time, data = sample_data)
  design_null <- model.matrix(~ replicate, data = sample_data)
  one_way <- FALSE
} else if (opt$model == "linear_time") {
  design_full <- model.matrix(~ time_numeric, data = sample_data)
  design_null <- model.matrix(~ 1, data = sample_data)
  one_way <- FALSE
} else {
  design_full <- model.matrix(~ replicate + time_numeric, data = sample_data)
  design_null <- model.matrix(~ replicate, data = sample_data)
  one_way <- FALSE
}
stopifnot(qr(design_full)$rank == ncol(design_full))

if (opt$threads > 1L) {
  bp <- BiocParallel::MulticoreParam(workers = opt$threads)
} else {
  bp <- BiocParallel::SerialParam()
}

set.seed(20260928)
d <- dmPrecision(d, design = design_full, one_way = one_way, BPPARAM = bp)
d <- dmFit(d, design = design_full, one_way = one_way, BPPARAM = bp)
d <- dmTest(d, design = design_null, one_way = one_way, BPPARAM = bp)

gene_results <- as.data.table(DRIMSeq::results(d, level = "gene"))
tx_results <- as.data.table(DRIMSeq::results(d, level = "feature"))
gene_results[is.na(pvalue), pvalue := 1]
tx_results[is.na(pvalue), pvalue := 1]

# stageR controls the overall false discovery rate across gene screening and
# transcript confirmation. The input p-values here are unadjusted.
p_screen <- gene_results$pvalue
names(p_screen) <- gene_results$gene_id
p_confirmation <- matrix(
  tx_results$pvalue,
  ncol = 1L,
  dimnames = list(tx_results$feature_id, "omnibus_time")
)
tx2gene <- as.data.frame(tx_results[, .(feature_id, gene_id)])
stage <- stageRTx(
  pScreen = p_screen,
  pConfirmation = p_confirmation,
  pScreenAdjusted = FALSE,
  tx2gene = tx2gene
)
stage <- stageWiseAdjustment(stage, method = "dtu", alpha = opt$alpha)
stage_table <- as.data.table(getAdjustedPValues(
  stage,
  onlySignificantGenes = FALSE,
  order = FALSE
))

# Calculate model-count usage among retained transcripts. This matches the
# composition tested by DRIMSeq but is not the only biological denominator.
retained_counts <- as.data.table(counts(d))
sample_cols <- samples$sample_id
gene_totals <- retained_counts[, lapply(.SD, sum), by = gene_id, .SDcols = sample_cols]
model_usage <- merge(retained_counts, gene_totals, by = "gene_id", suffixes = c("", "_gene"))
for (column in sample_cols) {
  denom <- paste0(column, "_gene")
  set(model_usage, j = column, value = fifelse(
    model_usage[[denom]] > 0,
    model_usage[[column]] / model_usage[[denom]],
    0
  ))
}
model_usage <- model_usage[, c("gene_id", "feature_id", sample_cols), with = FALSE]

for (minute in c(0L, 12L, 30L, 60L)) {
  cols <- samples[time_min == minute, sample_id]
  model_usage[, paste0("model_count_mean_usage_", minute) := rowMeans(.SD), .SDcols = cols]
}
for (minute in c(12L, 30L, 60L)) {
  model_usage[, paste0("model_count_delta_usage_", minute, "_vs_0") :=
          get(paste0("model_count_mean_usage_", minute)) - model_count_mean_usage_0]
}
model_delta_cols <- paste0("model_count_delta_usage_", c(12L, 30L, 60L), "_vs_0")
model_usage[, model_count_max_abs_delta_usage :=
              do.call(pmax, c(lapply(.SD, abs), list(na.rm = TRUE))),
            .SDcols = model_delta_cols]
model_usage[, (sample_cols) := NULL]

# Separately calculate count usage against every annotated transcript in the
# gene, including transcripts removed by the information filter.
all_count_totals <- count_table[, lapply(.SD, sum), by = gene_id, .SDcols = sample_cols]
all_count_usage <- merge(
  retained_counts[, .(gene_id, feature_id)],
  count_table,
  by = c("gene_id", "feature_id")
)
all_count_usage <- merge(all_count_usage, all_count_totals, by = "gene_id", suffixes = c("", "_gene"))
for (column in sample_cols) {
  denom <- paste0(column, "_gene")
  set(all_count_usage, j = column, value = fifelse(
    all_count_usage[[denom]] > 0,
    all_count_usage[[column]] / all_count_usage[[denom]],
    0
  ))
}
for (minute in c(0L, 12L, 30L, 60L)) {
  cols <- samples[time_min == minute, sample_id]
  all_count_usage[, paste0("all_count_mean_usage_", minute) := rowMeans(.SD), .SDcols = cols]
}
for (minute in c(12L, 30L, 60L)) {
  all_count_usage[, paste0("all_count_delta_usage_", minute, "_vs_0") :=
                    get(paste0("all_count_mean_usage_", minute)) - all_count_mean_usage_0]
}
all_count_delta_cols <- paste0("all_count_delta_usage_", c(12L, 30L, 60L), "_vs_0")
all_count_usage[, all_count_max_abs_delta_usage :=
                  do.call(pmax, c(lapply(.SD, abs), list(na.rm = TRUE))),
                .SDcols = all_count_delta_cols]
all_count_usage <- all_count_usage[, c(
  "gene_id", "feature_id", paste0("all_count_mean_usage_", c(0L, 12L, 30L, 60L)),
  all_count_delta_cols, "all_count_max_abs_delta_usage"
), with = FALSE]

# TPM fractions are closer to relative transcript molecule abundance and are
# therefore the primary descriptive effect size for the biological hypothesis.
tpm_list <- lapply(seq_len(nrow(samples)), function(i) {
  q <- as.data.table(read.delim(
    gzfile(file.path(project, samples$quant_file[i])),
    check.names = FALSE
  ))[, .(Name, TPM)]
  setnames(q, c("feature_id", samples$sample_id[i]))
  q
})
tpm_table <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), tpm_list)
tpm_table[is.na(tpm_table)] <- 0
tpm_table <- merge(
  annotation[, .(feature_id = transcript_id, gene_id)],
  tpm_table,
  by = "feature_id",
  all.y = TRUE
)
tpm_totals <- tpm_table[, lapply(.SD, sum), by = gene_id, .SDcols = sample_cols]
tpm_usage <- merge(retained_counts[, .(gene_id, feature_id)], tpm_table, by = c("gene_id", "feature_id"))
tpm_usage <- merge(tpm_usage, tpm_totals, by = "gene_id", suffixes = c("", "_gene"))
for (column in sample_cols) {
  denom <- paste0(column, "_gene")
  set(tpm_usage, j = column, value = fifelse(
    tpm_usage[[denom]] > 0,
    tpm_usage[[column]] / tpm_usage[[denom]],
    0
  ))
}
for (minute in c(0L, 12L, 30L, 60L)) {
  cols <- samples[time_min == minute, sample_id]
  tpm_usage[, paste0("tpm_mean_usage_", minute) := rowMeans(.SD), .SDcols = cols]
}
for (minute in c(12L, 30L, 60L)) {
  tpm_usage[, paste0("tpm_delta_usage_", minute, "_vs_0") :=
              get(paste0("tpm_mean_usage_", minute)) - tpm_mean_usage_0]
}
tpm_delta_cols <- paste0("tpm_delta_usage_", c(12L, 30L, 60L), "_vs_0")
tpm_usage[, tpm_max_abs_delta_usage :=
            do.call(pmax, c(lapply(.SD, abs), list(na.rm = TRUE))),
          .SDcols = tpm_delta_cols]
tpm_usage <- tpm_usage[, c(
  "gene_id", "feature_id", paste0("tpm_mean_usage_", c(0L, 12L, 30L, 60L)),
  tpm_delta_cols, "tpm_max_abs_delta_usage"
), with = FALSE]

tx_results <- merge(tx_results, model_usage, by = c("gene_id", "feature_id"), all.x = TRUE)
tx_results <- merge(tx_results, all_count_usage, by = c("gene_id", "feature_id"), all.x = TRUE)
tx_results <- merge(tx_results, tpm_usage, by = c("gene_id", "feature_id"), all.x = TRUE)
tx_results <- merge(
  tx_results,
  annotation,
  by.x = c("gene_id", "feature_id"),
  by.y = c("gene_id", "transcript_id"),
  all.x = TRUE
)

# stageR output format is retained as provided by the package and also merged
# into the transcript table after checking its identifier columns at runtime.
fwrite(gene_results, file.path(out_dir, "gene_omnibus.tsv"), sep = "\t")
fwrite(tx_results, file.path(out_dir, "transcript_omnibus_and_usage.tsv.gz"), sep = "\t")
fwrite(stage_table, file.path(out_dir, "stageR_adjusted.tsv"), sep = "\t")
fwrite(as.data.table(design_full, keep.rownames = "sample_id"), file.path(out_dir, "design_full.tsv"), sep = "\t")
fwrite(as.data.table(design_null, keep.rownames = "sample_id"), file.path(out_dir, "design_null.tsv"), sep = "\t")

summary <- data.table(
  filter = opt$filter,
  model = opt$model,
  input_scale = opt$input_scale,
  alpha = opt$alpha,
  retained_genes = uniqueN(tx_results$gene_id),
  retained_transcripts = nrow(tx_results),
  drimseq_gene_fdr_0_05 = sum(gene_results$adj_pvalue <= opt$alpha, na.rm = TRUE),
  drimseq_tx_fdr_0_05 = sum(tx_results$adj_pvalue <= opt$alpha, na.rm = TRUE)
)
fwrite(summary, file.path(out_dir, "run_summary.tsv"), sep = "\t")
saveRDS(d, file.path(out_dir, "drimseq_fit.rds"), compress = "xz")
print(summary)
print(head(stage_table))
