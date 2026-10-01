#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tximport)
  library(DESeq2)
})

project <- normalizePath(getwd(), mustWork = TRUE)
samples <- fread(file.path(project, "metadata", "samples.tsv"))
annotation <- as.data.table(read.delim(
  gzfile(file.path(project, "metadata", "gencode_v35_transcripts.tsv.gz")),
  check.names = FALSE
))
out_dir <- file.path(project, "results", "gene_expression_time")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

quant_files <- file.path(project, samples$quant_file)
names(quant_files) <- samples$sample_id
tx2gene <- unique(annotation[, .(transcript_id, gene_id)])

# Gene-level abundance and length correction are imported directly from the
# existing Salmon estimates. DESeq2 uses the tximport average-transcript-length
# offsets rather than treating transcript length as biological expression.
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = as.data.frame(tx2gene),
  dropInfReps = TRUE
)

coldata <- data.frame(
  sample_id = samples$sample_id,
  time = factor(samples$time_min, levels = c(0, 12, 30, 60)),
  row.names = samples$sample_id
)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ time)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
res <- as.data.table(results(dds), keep.rownames = "gene_id")
res <- merge(
  res,
  unique(annotation[, .(gene_id, gene_name, gene_type)]),
  by = "gene_id",
  all.x = TRUE
)

norm <- as.data.table(counts(dds, normalized = TRUE), keep.rownames = "gene_id")
for (minute in c(0L, 12L, 30L, 60L)) {
  cols <- samples[time_min == minute, sample_id]
  norm[, paste0("normalized_mean_", minute) := rowMeans(.SD), .SDcols = cols]
}
for (minute in c(12L, 30L, 60L)) {
  norm[, paste0("descriptive_log2fc_", minute, "_vs_0") :=
         log2((get(paste0("normalized_mean_", minute)) + 0.5) /
              (normalized_mean_0 + 0.5))]
}
norm[, max_abs_descriptive_log2fc := do.call(
  pmax,
  c(lapply(.SD, abs), list(na.rm = TRUE))
), .SDcols = paste0("descriptive_log2fc_", c(12L, 30L, 60L), "_vs_0")]

keep <- c(
  "gene_id", paste0("normalized_mean_", c(0L, 12L, 30L, 60L)),
  paste0("descriptive_log2fc_", c(12L, 30L, 60L), "_vs_0"),
  "max_abs_descriptive_log2fc"
)
res <- merge(res, norm[, ..keep], by = "gene_id", all.x = TRUE)
setorder(res, padj, pvalue)
fwrite(res, file.path(out_dir, "gene_time_omnibus_deseq2.tsv.gz"), sep = "\t")

summary <- data.table(
  tested_nonzero_genes = nrow(res),
  genes_with_padj = sum(!is.na(res$padj)),
  genes_padj_le_0_05 = sum(res$padj <= 0.05, na.rm = TRUE),
  model = "DESeq2 LRT: full ~ time; reduced ~ 1"
)
fwrite(summary, file.path(out_dir, "run_summary.tsv"), sep = "\t")
print(summary)
