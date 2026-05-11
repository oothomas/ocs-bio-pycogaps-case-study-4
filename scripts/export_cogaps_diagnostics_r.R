#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "."
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
default_preprocessed <- file.path(ROOT, "data", "processed", "preprocessed_cells_hvg3000.h5ad")
default_result_rds <- file.path(ROOT, "data", "results_r_k6_sparse_mt_t4_heavy", "cogaps_K6_seed2_iter2000.rds")
default_out_dir <- file.path(ROOT, "data", "results_r_k6_sparse_mt_t4_heavy")

option_list <- list(
  make_option("--preprocessed-h5ad", type = "character", dest = "preprocessed_h5ad", default = default_preprocessed, help = "Processed cells x genes AnnData (.h5ad) [default %default]"),
  make_option("--result-rds", type = "character", dest = "result_rds", default = default_result_rds, help = "Saved CoGAPS R result (.rds) [default %default]"),
  make_option("--outdir", type = "character", default = default_out_dir, help = "Output directory [default %default]"),
  make_option("--top-n", type = "integer", dest = "top_n", default = 15, help = "Top features used in stability summaries [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

PREPROCESSED_H5AD <- opt$preprocessed_h5ad
RESULT_RDS <- opt$result_rds
OUT_DIR <- opt$outdir
TOP_N <- opt$top_n

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

tag <- sub("^cogaps_(.*)\\.rds$", "\\1", basename(RESULT_RDS))
if (identical(tag, basename(RESULT_RDS))) {
  tag <- tools::file_path_sans_ext(basename(RESULT_RDS))
}

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

write_gz_tsv <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

safe_cor <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2) {
    return(NA_real_)
  }
  if (sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y))
}

top_index_overlap <- function(current, final, n = 15) {
  n <- min(n, length(current), length(final))
  current_top <- names(sort(current, decreasing = TRUE))[seq_len(n)]
  final_top <- names(sort(final, decreasing = TRUE))[seq_len(n)]
  length(intersect(current_top, final_top)) / n
}

param_to_list <- function(params) {
  out <- list()
  for (nm in slotNames(params)) {
    value <- slot(params, nm)
    if (is.null(value)) {
      out[[nm]] <- NULL
    } else if (length(dim(value)) > 0) {
      out[[nm]] <- list(dim = dim(value), class = class(value))
    } else if (is.atomic(value) && length(value) <= 50) {
      out[[nm]] <- unname(value)
    } else {
      out[[nm]] <- list(class = class(value), length = length(value))
    }
  }
  out
}

build_trace_df <- function(md) {
  n <- max(length(md$chisq %||% numeric()), length(md$atomsA %||% numeric()), length(md$atomsP %||% numeric()))
  if (!is.finite(n) || n == 0) {
    return(data.frame())
  }
  fill <- function(x) {
    x <- x %||% numeric()
    c(as.numeric(x), rep(NA_real_, n - length(x)))
  }
  out <- data.frame(
    trace_index = seq_len(n),
    relative_progress = seq_len(n) / n,
    chisq = fill(md$chisq),
    atomsA = fill(md$atomsA),
    atomsP = fill(md$atomsP)
  )
  out$delta_chisq <- c(NA_real_, diff(out$chisq))
  out$delta_atomsA <- c(NA_real_, diff(out$atomsA))
  out$delta_atomsP <- c(NA_real_, diff(out$atomsP))
  out
}

summarize_snapshots <- function(snapshot_list, final_mat, phase, matrix_name, top_n) {
  if (!length(snapshot_list)) {
    return(data.frame())
  }
  final_vec <- as.numeric(final_mat)
  final_norm <- sqrt(sum(final_vec^2))
  do.call(
    rbind,
    lapply(seq_along(snapshot_list), function(i) {
      snap <- snapshot_list[[i]]
      if (!is.null(colnames(snap))) {
        colnames(snap) <- canonicalize_pattern_names(colnames(snap))
      }
      snap_vec <- as.numeric(snap)
      overlap <- NA_real_
      if (identical(matrix_name, "A")) {
        per_pattern_overlap <- vapply(
          colnames(final_mat),
          function(pattern) top_index_overlap(snap[, pattern], final_mat[, pattern], n = top_n),
          numeric(1)
        )
        overlap <- mean(per_pattern_overlap)
      }
      data.frame(
        phase = phase,
        snapshot_label = names(snapshot_list)[i] %||% as.character(i),
        matrix = matrix_name,
        frobenius_to_final = sqrt(sum((snap_vec - final_vec)^2)),
        relative_frobenius_to_final = sqrt(sum((snap_vec - final_vec)^2)) / final_norm,
        correlation_to_final = safe_cor(snap_vec, final_vec),
        nnz = sum(snap_vec > 0),
        mean_top_overlap = overlap,
        row.names = NULL
      )
    })
  )
}

summarize_uncertainty <- function(A, P, A_sd, P_sd, top_n) {
  patterns <- colnames(A)
  do.call(
    rbind,
    lapply(patterns, function(pattern) {
      ordered_genes <- order(A[, pattern], decreasing = TRUE)
      top_idx <- ordered_genes[seq_len(min(top_n, length(ordered_genes)))]
      weight_top <- A[top_idx, pattern]
      sd_top <- A_sd[top_idx, pattern]
      cv_top <- sd_top / pmax(abs(weight_top), 1e-8)
      data.frame(
        pattern = pattern,
        mean_loading_sd = mean(A_sd[, pattern]),
        median_loading_sd = median(A_sd[, pattern]),
        max_loading_sd = max(A_sd[, pattern]),
        mean_factor_sd = mean(P_sd[, pattern]),
        median_factor_sd = median(P_sd[, pattern]),
        max_factor_sd = max(P_sd[, pattern]),
        top_gene_mean_weight = mean(weight_top),
        top_gene_mean_loading_sd = mean(sd_top),
        top_gene_median_loading_cv = median(cv_top),
        row.names = NULL
      )
    })
  )
}

sce <- zellkonverter::readH5AD(PREPROCESSED_H5AD, reader = "R")
result <- readRDS(RESULT_RDS)
md <- slot(result, "metadata")

A <- as.matrix(getFeatureLoadings(result))
P <- as.matrix(getSampleFactors(result))
A_sd <- as.matrix(slot(result, "loadingStdDev"))
P_sd <- as.matrix(slot(result, "factorStdDev"))

A <- ensure_pattern_matrix_orientation(A, nrow(sce), rownames(sce))
P <- ensure_pattern_matrix_orientation(P, ncol(sce), colnames(sce))
A_sd <- ensure_pattern_matrix_orientation(A_sd, nrow(sce), rownames(sce))
P_sd <- ensure_pattern_matrix_orientation(P_sd, ncol(sce), colnames(sce))

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
colnames(A_sd) <- colnames(A)
colnames(P_sd) <- colnames(P)

trace_df <- build_trace_df(md)
snapshot_summary <- rbind(
  summarize_snapshots(md$equilibrationSnapshotsA %||% list(), A, "equilibration", "A", TOP_N),
  summarize_snapshots(md$equilibrationSnapshotsP %||% list(), P, "equilibration", "P", TOP_N),
  summarize_snapshots(md$samplingSnapshotsA %||% list(), A, "sampling", "A", TOP_N),
  summarize_snapshots(md$samplingSnapshotsP %||% list(), P, "sampling", "P", TOP_N)
)
uncertainty_summary <- summarize_uncertainty(A, P, A_sd, P_sd, top_n = TOP_N)

if (!is.null(md$pumpStat) && all(dim(md$pumpStat) > 0)) {
  pump_df <- data.frame(gene = rownames(A), md$pumpStat, row.names = NULL, check.names = FALSE)
  colnames(pump_df)[-1] <- colnames(A)
  write_gz_tsv(pump_df, file.path(OUT_DIR, sprintf("cogaps_%s.pump_stat.tsv.gz", tag)))
}

if (!is.null(md$meanPatternAssignment) && all(dim(md$meanPatternAssignment) > 0)) {
  assign_df <- data.frame(gene = rownames(A), md$meanPatternAssignment, row.names = NULL, check.names = FALSE)
  colnames(assign_df)[-1] <- colnames(A)
  write_gz_tsv(assign_df, file.path(OUT_DIR, sprintf("cogaps_%s.mean_pattern_assignment.tsv.gz", tag)))
}

write.csv(trace_df, file.path(OUT_DIR, sprintf("cogaps_%s.trace.csv", tag)), row.names = FALSE)
write.csv(snapshot_summary, file.path(OUT_DIR, sprintf("cogaps_%s.snapshot_summary.csv", tag)), row.names = FALSE)
write.csv(uncertainty_summary, file.path(OUT_DIR, sprintf("cogaps_%s.pattern_uncertainty_summary.csv", tag)), row.names = FALSE)

diagnostics_summary <- list(
  status = "ok",
  language = "R",
  package = "CoGAPS",
  package_version = as.character(md$version %||% utils::packageVersion("CoGAPS")),
  result_rds = basename(RESULT_RDS),
  dimensions = list(genes = nrow(A), samples = nrow(P), patterns = ncol(A)),
  params = param_to_list(md$params),
  metadata_summary = list(
    meanChiSq = md$meanChiSq %||% NA_real_,
    totalUpdates = md$totalUpdates %||% NA_real_,
    totalRunningTime = md$totalRunningTime %||% NA_real_,
    averageQueueLengthA = md$averageQueueLengthA %||% NA_real_,
    averageQueueLengthP = md$averageQueueLengthP %||% NA_real_,
    chisq_points = length(md$chisq %||% numeric()),
    atomsA_points = length(md$atomsA %||% numeric()),
    atomsP_points = length(md$atomsP %||% numeric()),
    equilibrationSnapshots = length(md$equilibrationSnapshotsA %||% list()),
    samplingSnapshots = length(md$samplingSnapshotsA %||% list()),
    pumpStatDim = dim(md$pumpStat %||% matrix(numeric(), nrow = 0, ncol = 0)),
    meanPatternAssignmentDim = dim(md$meanPatternAssignment %||% matrix(numeric(), nrow = 0, ncol = 0))
  ),
  trace_summary = if (nrow(trace_df)) list(
    chisq_start = trace_df$chisq[1],
    chisq_end = trace_df$chisq[nrow(trace_df)],
    chisq_min = min(trace_df$chisq, na.rm = TRUE),
    chisq_max = max(trace_df$chisq, na.rm = TRUE),
    chisq_num_increases = sum(diff(trace_df$chisq) > 0, na.rm = TRUE),
    chisq_num_decreases = sum(diff(trace_df$chisq) < 0, na.rm = TRUE),
    atomsA_start = trace_df$atomsA[1],
    atomsA_end = trace_df$atomsA[nrow(trace_df)],
    atomsP_start = trace_df$atomsP[1],
    atomsP_end = trace_df$atomsP[nrow(trace_df)]
  ) else NULL
)

write_json(
  diagnostics_summary,
  path = file.path(OUT_DIR, sprintf("cogaps_%s.diagnostics.json", tag)),
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null"
)

cat("Wrote CoGAPS diagnostics to:", OUT_DIR, "\n")
