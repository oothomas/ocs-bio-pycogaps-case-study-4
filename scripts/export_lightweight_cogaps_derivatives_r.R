#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "."
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
PREPROCESSED_H5AD <- file.path(ROOT, "data", "processed", "preprocessed_cells_hvg3000.h5ad")
CHOSEN_RESULT_RDS <- file.path(ROOT, "data", "results_r", "cogaps_K7_seed2_iter2000.rds")
OUT_DIR <- file.path(ROOT, "data", "results_r")
STIM_LABEL <- "stim"

pattern_columns <- function(x) {
  pats <- grep("^Pattern_?[0-9]+$", x, value = TRUE)
  pats[order(as.integer(sub("^Pattern_?", "", pats)))]
}

canonicalize_pattern_names <- function(x) {
  sub("^Pattern_?([0-9]+)$", "Pattern\\1", x)
}

ensure_pattern_matrix_orientation <- function(mat, expected_rows, expected_names) {
  if (nrow(mat) != expected_rows && ncol(mat) == expected_rows) {
    mat <- t(mat)
  }
  if (nrow(mat) != expected_rows) {
    stop(sprintf("Unexpected matrix shape. Expected %d rows, got %d.", expected_rows, nrow(mat)))
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- expected_names
  }
  mat
}

top_genes_from_matrix <- function(matrix, pattern, n = 50) {
  ordered <- sort(matrix[, pattern], decreasing = TRUE)
  ordered <- head(ordered, n)
  data.frame(
    pattern = pattern,
    rank = seq_along(ordered),
    gene = names(ordered),
    weight = unname(ordered),
    row.names = NULL
  )
}

write_gz_tsv <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

coerce_aggregate_stats <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    return(as.matrix(x))
  }
  do.call(rbind, x)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sce <- zellkonverter::readH5AD(PREPROCESSED_H5AD, reader = "R")
result <- readRDS(CHOSEN_RESULT_RDS)

A <- as.matrix(getFeatureLoadings(result))
P <- as.matrix(getSampleFactors(result))

A <- ensure_pattern_matrix_orientation(A, nrow(sce), rownames(sce))
P <- ensure_pattern_matrix_orientation(P, ncol(sce), colnames(sce))

if (is.null(colnames(A))) {
  colnames(A) <- sprintf("Pattern%d", seq_len(ncol(A)))
} else {
  colnames(A) <- canonicalize_pattern_names(colnames(A))
}
if (is.null(colnames(P))) {
  colnames(P) <- sprintf("Pattern%d", seq_len(ncol(P)))
} else {
  colnames(P) <- canonicalize_pattern_names(colnames(P))
}

pattern_names <- pattern_columns(colnames(A))

cell_meta <- as.data.frame(colData(sce))
cell_meta$cell_barcode <- colnames(sce)

if ("X_umap" %in% reducedDimNames(sce)) {
  umap <- reducedDim(sce, "X_umap")
  cell_meta$umap_1 <- umap[, 1]
  cell_meta$umap_2 <- umap[, 2]
}

if ("X_pca" %in% reducedDimNames(sce)) {
  pca <- reducedDim(sce, "X_pca")
  cell_meta$pca_1 <- pca[, 1]
  cell_meta$pca_2 <- pca[, 2]
}

P_df <- as.data.frame(P)
P_df$cell_barcode <- rownames(P)
cell_activities_with_metadata <- merge(cell_meta, P_df, by = "cell_barcode", sort = FALSE)

stim_indicator <- as.integer(as.character(cell_meta$condition) == STIM_LABEL)
pattern_corrs <- data.frame(
  pattern = pattern_names,
  corr_with_stim = vapply(pattern_names, function(pattern) cor(P[, pattern], stim_indicator), numeric(1)),
  row.names = NULL
)
pattern_corrs <- pattern_corrs[order(pattern_corrs$corr_with_stim, decreasing = TRUE), ]

top_genes <- do.call(rbind, lapply(pattern_names, function(pattern) top_genes_from_matrix(A, pattern, n = 50)))

activity_by_celltype_condition <- do.call(
  rbind,
  lapply(pattern_names, function(pattern) {
    ag <- aggregate(
      cell_activities_with_metadata[[pattern]],
      by = list(
        cell_type = cell_activities_with_metadata$cell_type,
        condition = cell_activities_with_metadata$condition
      ),
      FUN = function(x) c(mean = mean(x), median = median(x), n_cells = length(x))
    )
    stats <- coerce_aggregate_stats(ag$x)
    data.frame(
      pattern = pattern,
      cell_type = ag$cell_type,
      condition = ag$condition,
      mean_activity = stats[, "mean"],
      median_activity = stats[, "median"],
      n_cells = stats[, "n_cells"],
      row.names = NULL
    )
  })
)

activity_by_replicate_condition <- do.call(
  rbind,
  lapply(pattern_names, function(pattern) {
    ag <- aggregate(
      cell_activities_with_metadata[[pattern]],
      by = list(
        replicate = cell_activities_with_metadata$replicate,
        condition = cell_activities_with_metadata$condition
      ),
      FUN = function(x) c(mean = mean(x), median = median(x), n_cells = length(x))
    )
    stats <- coerce_aggregate_stats(ag$x)
    data.frame(
      pattern = pattern,
      replicate = ag$replicate,
      condition = ag$condition,
      mean_activity = stats[, "mean"],
      median_activity = stats[, "median"],
      n_cells = stats[, "n_cells"],
      row.names = NULL
    )
  })
)

pattern_summary <- do.call(
  rbind,
  lapply(pattern_names, function(pattern) {
    top8 <- head(top_genes$gene[top_genes$pattern == pattern], 8)
    top15 <- head(top_genes$gene[top_genes$pattern == pattern], 15)
    celltype_means <- tapply(cell_activities_with_metadata[[pattern]], cell_activities_with_metadata$cell_type, mean)
    condition_means <- tapply(cell_activities_with_metadata[[pattern]], cell_activities_with_metadata$condition, mean)
    data.frame(
      pattern = pattern,
      corr_with_stim = pattern_corrs$corr_with_stim[pattern_corrs$pattern == pattern],
      dominant_cell_type = names(sort(celltype_means, decreasing = TRUE))[1],
      dominant_cell_type_mean_activity = unname(sort(celltype_means, decreasing = TRUE)[1]),
      mean_activity_ctrl = unname(condition_means["ctrl"]),
      mean_activity_stim = unname(condition_means["stim"]),
      top_genes_top8 = paste(top8, collapse = ", "),
      top_genes_top15 = paste(top15, collapse = ", "),
      row.names = NULL
    )
  })
)

A_df <- data.frame(gene = rownames(A), A, row.names = NULL, check.names = FALSE)
P_df_out <- data.frame(cell_barcode = rownames(P), P, row.names = NULL, check.names = FALSE)

write_gz_tsv(A_df, file.path(OUT_DIR, "pattern_gene_weights.tsv.gz"))
write_gz_tsv(P_df_out, file.path(OUT_DIR, "pattern_cell_activities.tsv.gz"))
write.csv(cell_activities_with_metadata, file = gzfile(file.path(OUT_DIR, "pattern_cell_activities_with_metadata.csv.gz"), "wt"), row.names = FALSE)
write.csv(top_genes, file.path(OUT_DIR, "pattern_top_genes.csv"), row.names = FALSE)
write.csv(pattern_corrs, file.path(OUT_DIR, "pattern_correlations.csv"), row.names = FALSE)
write.csv(activity_by_celltype_condition, file.path(OUT_DIR, "pattern_activity_by_celltype_condition.csv"), row.names = FALSE)
write.csv(activity_by_replicate_condition, file.path(OUT_DIR, "pattern_activity_by_replicate_condition.csv"), row.names = FALSE)
write.csv(pattern_summary, file.path(OUT_DIR, "pattern_summary.csv"), row.names = FALSE)

cat("Wrote lightweight R CoGAPS derivatives to:", OUT_DIR, "\n")
