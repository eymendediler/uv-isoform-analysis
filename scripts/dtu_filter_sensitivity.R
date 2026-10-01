#!/usr/bin/env Rscript

# Compare reasonable low-information filters without performing a DTU test.
# The purpose is to inspect data retention before jointly choosing a primary
# analysis filter. No sample is removed and no significance cutoff is applied.

suppressPackageStartupMessages({
  library(DRIMSeq)
  library(data.table)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
project <- normalizePath(file.path(dirname(file_arg), ".."))

sample_file <- file.path(project, "metadata", "samples.tsv")
annotation_file <- file.path(project, "metadata", "gencode_v35_transcripts.tsv.gz")
out_dir <- file.path(project, "results", "dtu_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- fread(sample_file)
stopifnot(nrow(samples) == 12L, all(table(samples$time_min) == 3L))
annotation <- as.data.table(read.delim(gzfile(annotation_file), check.names = FALSE))

count_list <- lapply(seq_len(nrow(samples)), function(i) {
  q <- as.data.table(read.delim(
    gzfile(file.path(project, samples$quant_file[i])),
    check.names = FALSE
  ))[, .(Name, NumReads)]
  setnames(q, c("feature_id", samples$sample_id[i]))
  q
})
counts <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), count_list)
counts[is.na(counts)] <- 0
for (column in samples$sample_id) set(counts, j = column, value = as.integer(round(counts[[column]])))

counts <- merge(
  annotation[, .(feature_id = transcript_id, gene_id)],
  counts,
  by = "feature_id",
  all.y = TRUE
)
stopifnot(!anyNA(counts$gene_id), nrow(counts) == 229580L)
setcolorder(counts, c("gene_id", "feature_id", samples$sample_id))

sample_data <- data.frame(
  sample_id = samples$sample_id,
  time = factor(samples$time_min, levels = c(0, 12, 30, 60)),
  replicate = factor(samples$replicate),
  row.names = samples$sample_id,
  check.names = FALSE
)
d0 <- dmDSdata(counts = as.data.frame(counts), samples = sample_data)

grid <- data.table(
  filter = c("loose", "medium", "strict"),
  min_samps_gene_expr = c(3L, 3L, 3L),
  min_gene_expr = c(10L, 20L, 50L),
  min_samps_feature_expr = c(2L, 3L, 3L),
  min_feature_expr = c(5L, 10L, 20L),
  min_samps_feature_prop = c(2L, 3L, 3L),
  min_feature_prop = c(0.01, 0.05, 0.10)
)

retention <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i]
  filtered <- dmFilter(
    d0,
    min_samps_gene_expr = g$min_samps_gene_expr,
    min_gene_expr = g$min_gene_expr,
    min_samps_feature_expr = g$min_samps_feature_expr,
    min_feature_expr = g$min_feature_expr,
    min_samps_feature_prop = g$min_samps_feature_prop,
    min_feature_prop = g$min_feature_prop
  )
  retained <- counts(filtered)
  data.table(
    filter = g$filter,
    retained_genes = length(unique(retained$gene_id)),
    retained_transcripts = nrow(retained),
    genes_with_multiple_retained_transcripts = sum(table(retained$gene_id) >= 2L)
  )
}))

report <- merge(grid, retention, by = "filter")
report[, original_genes := uniqueN(counts$gene_id)]
report[, original_transcripts := nrow(counts)]
report[, gene_retention_fraction := retained_genes / original_genes]
report[, transcript_retention_fraction := retained_transcripts / original_transcripts]
setcolorder(report, c(
  "filter", "min_samps_gene_expr", "min_gene_expr",
  "min_samps_feature_expr", "min_feature_expr",
  "min_samps_feature_prop", "min_feature_prop",
  "original_genes", "retained_genes", "gene_retention_fraction",
  "original_transcripts", "retained_transcripts", "transcript_retention_fraction",
  "genes_with_multiple_retained_transcripts"
))
fwrite(report, file.path(out_dir, "filter_retention.tsv"), sep = "\t")
print(report)
