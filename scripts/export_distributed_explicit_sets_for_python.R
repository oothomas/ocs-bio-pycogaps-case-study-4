#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(zellkonverter)
})

option_list <- list(
  make_option("--explicit-sets-rds", type = "character", dest = "explicit_sets_rds", help = "RDS file containing the R explicit sets"),
  make_option("--preprocessed-h5ad", type = "character", dest = "preprocessed_h5ad", default = NULL, help = "Optional processed cells x genes H5AD used to attach cell IDs"),
  make_option("--out-json", type = "character", dest = "out_json", help = "Output JSON path"),
  make_option("--subset-summary-csv", type = "character", dest = "subset_summary_csv", default = NULL, help = "Optional output CSV with set sizes"),
  make_option("--label", type = "character", dest = "label", default = "distributed_single_cell", help = "Free-text label stored in JSON metadata")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$explicit_sets_rds) || is.null(opt$out_json)) {
  stop("Both --explicit-sets-rds and --out-json are required.")
}

explicit_sets <- readRDS(opt$explicit_sets_rds)
if (!is.list(explicit_sets) || length(explicit_sets) == 0) {
  stop("The explicit sets RDS did not contain a non-empty list.")
}

cell_ids <- NULL
if (!is.null(opt$preprocessed_h5ad) && nzchar(opt$preprocessed_h5ad)) {
  sce <- zellkonverter::readH5AD(opt$preprocessed_h5ad, reader = "R")
  cell_ids <- colnames(sce)
}

sets_payload <- lapply(seq_along(explicit_sets), function(i) {
  one_based <- as.integer(explicit_sets[[i]])
  zero_based <- one_based - 1L
  if (any(zero_based < 0)) {
    stop("Encountered negative zero-based indices during export.")
  }

  out <- list(
    set_id = sprintf("set%02d", i),
    one_based_indices = unname(one_based),
    zero_based_indices = unname(zero_based)
  )
  if (!is.null(cell_ids)) {
    out$cell_ids <- unname(cell_ids[one_based])
  }
  out
})

payload <- list(
  status = "ok",
  label = opt$label,
  indexing = "zero-based",
  n_sets = length(explicit_sets),
  source_rds = normalizePath(opt$explicit_sets_rds, winslash = "/", mustWork = TRUE),
  preprocessed_h5ad = if (!is.null(opt$preprocessed_h5ad) && nzchar(opt$preprocessed_h5ad)) {
    normalizePath(opt$preprocessed_h5ad, winslash = "/", mustWork = TRUE)
  } else {
    NULL
  },
  exported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  sets = sets_payload
)

dir.create(dirname(opt$out_json), recursive = TRUE, showWarnings = FALSE)
write_json(payload, path = opt$out_json, auto_unbox = TRUE, pretty = TRUE)

if (!is.null(opt$subset_summary_csv) && nzchar(opt$subset_summary_csv)) {
  rows <- lapply(sets_payload, function(entry) {
    data.frame(
      set_id = entry$set_id,
      n_cells = length(entry$zero_based_indices),
      first_index = if (length(entry$zero_based_indices) > 0) entry$zero_based_indices[[1]] else NA_integer_,
      last_index = if (length(entry$zero_based_indices) > 0) tail(entry$zero_based_indices, 1) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, rows)
  dir.create(dirname(opt$subset_summary_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, opt$subset_summary_csv, row.names = FALSE)
}
