#!/usr/bin/env Rscript

meta <- read.csv(
  gzfile("data/results_r_sparse_instrumented/pattern_cell_activities_with_metadata.csv.gz"),
  check.names = FALSE
)
meta <- meta[, c("cell_barcode", "cell_type", "condition", "replicate")]

export_activity <- function(rds, out_csv) {
  x <- readRDS(rds)
  P <- as.data.frame(x@sampleFactors, check.names = FALSE)
  colnames(P) <- gsub("Pattern_", "Pattern", colnames(P), fixed = TRUE)
  P$cell_barcode <- rownames(x@sampleFactors)
  df <- merge(meta, P, by = "cell_barcode", all.y = TRUE)
  pats <- grep("^Pattern", names(df), value = TRUE)

  rows <- list()
  i <- 1L
  for (pat in pats) {
    tmp <- aggregate(
      df[[pat]],
      by = list(pattern = rep(pat, nrow(df)), cell_type = df$cell_type, condition = df$condition),
      FUN = function(v) c(mean = mean(v, na.rm = TRUE), median = median(v, na.rm = TRUE), n = sum(!is.na(v)))
    )
    vals <- if (is.matrix(tmp$x)) tmp$x else do.call(rbind, tmp$x)
    rows[[i]] <- data.frame(
      pattern = tmp$pattern,
      cell_type = tmp$cell_type,
      condition = tmp$condition,
      mean_activity = vals[, "mean"],
      median_activity = vals[, "median"],
      n_cells = vals[, "n"],
      check.names = FALSE
    )
    i <- i + 1L
  }

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
}

export_activity(
  "data/results_r_distributed_sparse_on_full/cogaps_K7_seed2_iter2000_single_cell.rds",
  "data/diagnostics_model_comparison/tables/r_distributed_sparse_on_activity_by_celltype_condition.csv"
)
export_activity(
  "data/results_r_distributed_sparse_off_full/cogaps_K7_seed2_iter2000_single_cell.rds",
  "data/diagnostics_model_comparison/tables/r_distributed_sparse_off_activity_by_celltype_condition.csv"
)
