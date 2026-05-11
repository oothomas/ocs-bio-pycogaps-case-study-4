#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
})

pattern_columns <- function(x) {
  pats <- grep("^Pattern_?[0-9]+$", x, value = TRUE)
  pats[order(as.integer(sub("^Pattern_?", "", pats)))]
}

canonicalize_pattern_names <- function(x) {
  sub("^Pattern_?([0-9]+)$", "Pattern\\1", x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

ensure_pattern_matrix_orientation <- function(mat, expected_rows, expected_names) {
  if (is.null(dim(mat))) {
    stop("Expected a 2D matrix.")
  }

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

write_metrics <- function(path, metrics) {
  write_json(metrics, path = path, auto_unbox = TRUE, pretty = TRUE)
}

option_list <- list(
  make_option("--preprocessed-h5ad", type = "character", dest = "preprocessed_h5ad", help = "Processed cells x genes AnnData (.h5ad)"),
  make_option("--outdir", type = "character", default = "data/results_r", help = "Output directory [default %default]"),
  make_option("--k", type = "integer", help = "Number of patterns"),
  make_option("--seed", type = "integer", help = "Random seed"),
  make_option("--n-iter", type = "integer", dest = "n_iter", help = "Number of iterations"),
  make_option("--top-genes", type = "integer", dest = "top_genes", default = 50, help = "Top genes used for summary metrics [default %default]"),
  make_option("--stim-label", type = "character", dest = "stim_label", default = "stim", help = "Stimulated condition label [default %default]"),
  make_option("--cogaps-threads", type = "integer", dest = "cogaps_threads", default = 1, help = "Threads passed to CoGAPS [default %default]"),
  make_option("--async-updates", action = "store_true", dest = "asynchronous_updates", default = TRUE, help = "Enable CoGAPS asynchronous updates [default]"),
  make_option("--sync-updates", action = "store_false", dest = "asynchronous_updates", help = "Disable CoGAPS asynchronous updates"),
  make_option("--output-frequency", type = "integer", dest = "output_frequency", default = 1000, help = "Trace output frequency passed to CoGAPS [default %default]"),
  make_option("--checkpoint-interval", type = "integer", dest = "checkpoint_interval", default = 0, help = "Iterations between checkpoint writes [default %default]"),
  make_option("--checkpoint-outfile", type = "character", dest = "checkpoint_outfile", default = NULL, help = "Checkpoint output file path [default auto in outdir when checkpointing is enabled]"),
  make_option("--n-snapshots", type = "integer", dest = "n_snapshots", default = 0, help = "Snapshots to capture per enabled phase [default %default]"),
  make_option("--snapshot-phase", type = "character", dest = "snapshot_phase", default = "sampling", help = "Snapshot phase: sampling, equilibration, or all [default %default]"),
  make_option("--take-pump-samples", action = "store_true", dest = "take_pump_samples", default = FALSE, help = "Capture pump statistics and mean pattern assignment diagnostics"),
  make_option("--quiet", action = "store_false", dest = "messages", default = TRUE, help = "Suppress CoGAPS status messages"),
  make_option("--use-sparse-opt", action = "store_true", dest = "use_sparse_opt", default = TRUE, help = "Enable sparseOptimization [default]"),
  make_option("--no-sparse-opt", action = "store_false", dest = "use_sparse_opt", help = "Disable sparseOptimization"),
  make_option("--force-rerun", action = "store_true", default = FALSE, help = "Ignore successful cached outputs and rerun")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$preprocessed_h5ad) || is.null(opt$k) || is.null(opt$seed) || is.null(opt$n_iter)) {
  stop("Missing required arguments. Use --help for details.")
}

outdir <- normalizePath(opt$outdir, winslash = "/", mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tag <- sprintf("K%d_seed%d_iter%d", opt$k, opt$seed, opt$n_iter)
result_path <- file.path(outdir, sprintf("cogaps_%s.rds", tag))
metrics_path <- file.path(outdir, sprintf("cogaps_%s.metrics.json", tag))
checkpoint_path <- opt$checkpoint_outfile
if (is.null(checkpoint_path) || identical(checkpoint_path, "")) {
  checkpoint_path <- file.path(outdir, sprintf("cogaps_%s.checkpoint.out", tag))
}

if (file.exists(result_path) && file.exists(metrics_path) && !isTRUE(opt$force_rerun)) {
  prior_metrics <- tryCatch(fromJSON(metrics_path), error = function(e) NULL)
  if (!is.null(prior_metrics) && identical(prior_metrics$status, "ok")) {
    message(sprintf("[CACHE] status=ok; skipping %s", tag))
    quit(save = "no", status = 0)
  }
}

status <- "ok"
err_txt <- NULL
t0 <- Sys.time()

tryCatch({
  sce <- zellkonverter::readH5AD(opt$preprocessed_h5ad, reader = "R")
  data_matrix <- SummarizedExperiment::assay(sce, "X")
  data_matrix <- as.matrix(data_matrix)
  openmp_support <- compiledWithOpenMPSupport()

  if (opt$cogaps_threads > 1 && !openmp_support) {
    stop(
      paste(
        "This CoGAPS build does not report OpenMP support, so requesting",
        sprintf("nThreads=%d", opt$cogaps_threads),
        "would silently fall back to a single-threaded run.",
        "Rebuild the runtime with OpenMP-enabled CoGAPS before benchmarking multithreading."
      )
    )
  }

  params <- CogapsParams(
    nPatterns = opt$k,
    nIterations = opt$n_iter,
    seed = opt$seed,
    sparseOptimization = isTRUE(opt$use_sparse_opt)
  )
  params@takePumpSamples <- isTRUE(opt$take_pump_samples)

  set.seed(opt$seed)
  result <- CoGAPS(
    data_matrix,
    params,
    nThreads = opt$cogaps_threads,
    messages = isTRUE(opt$messages),
    outputFrequency = opt$output_frequency,
    checkpointOutFile = checkpoint_path,
    checkpointInterval = opt$checkpoint_interval,
    nSnapshots = opt$n_snapshots,
    snapshotPhase = opt$snapshot_phase,
    asynchronousUpdates = isTRUE(opt$asynchronous_updates)
  )
  saveRDS(result, result_path)

  feature_loadings <- as.matrix(getFeatureLoadings(result))
  sample_factors <- as.matrix(getSampleFactors(result))

  feature_loadings <- ensure_pattern_matrix_orientation(
    mat = feature_loadings,
    expected_rows = nrow(sce),
    expected_names = rownames(sce)
  )
  sample_factors <- ensure_pattern_matrix_orientation(
    mat = sample_factors,
    expected_rows = ncol(sce),
    expected_names = colnames(sce)
  )

  if (is.null(colnames(feature_loadings))) {
    colnames(feature_loadings) <- sprintf("Pattern%d", seq_len(ncol(feature_loadings)))
  } else {
    colnames(feature_loadings) <- canonicalize_pattern_names(colnames(feature_loadings))
  }
  if (is.null(colnames(sample_factors))) {
    colnames(sample_factors) <- sprintf("Pattern%d", seq_len(ncol(sample_factors)))
  } else {
    colnames(sample_factors) <- canonicalize_pattern_names(colnames(sample_factors))
  }

  pat_names <- pattern_columns(colnames(sample_factors))
  if (!"condition" %in% colnames(colData(sce))) {
    stop("The processed H5AD does not contain a 'condition' column in cell metadata.")
  }

  stim_indicator <- as.integer(as.character(colData(sce)$condition) == opt$stim_label)
  pattern_corrs <- vapply(
    pat_names,
    function(pattern) cor(sample_factors[, pattern], stim_indicator),
    numeric(1)
  )
  ifn_pattern <- names(which.max(pattern_corrs))
  ifn_corr <- unname(pattern_corrs[[ifn_pattern]])

  top_genes <- rownames(feature_loadings)[order(feature_loadings[, ifn_pattern], decreasing = TRUE)][seq_len(min(opt$top_genes, nrow(feature_loadings)))]

  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  md <- slot(result, "metadata")

  metrics <- list(
    status = "ok",
    language = "R",
    package = "CoGAPS",
    result_path = basename(result_path),
    nPatterns = opt$k,
    nIterations = opt$n_iter,
    seed = opt$seed,
    sparseOptimization = isTRUE(opt$use_sparse_opt),
    nThreads = opt$cogaps_threads,
    asynchronousUpdates = isTRUE(opt$asynchronous_updates),
    openmpSupport = openmp_support,
    outputFrequency = opt$output_frequency,
    checkpointInterval = opt$checkpoint_interval,
    checkpointOutFile = basename(checkpoint_path),
    nSnapshots = opt$n_snapshots,
    snapshotPhase = opt$snapshot_phase,
    takePumpSamples = isTRUE(opt$take_pump_samples),
    ifn_pattern = ifn_pattern,
    ifn_corr = ifn_corr,
    top_genes = unname(top_genes),
    runtime_sec = runtime_sec,
    totalRunningTime = md$totalRunningTime %||% NA_real_,
    totalUpdates = md$totalUpdates %||% NA_real_,
    meanChiSq = md$meanChiSq %||% NA_real_,
    chisq_trace_length = length(md$chisq %||% numeric()),
    atomsA_trace_length = length(md$atomsA %||% numeric()),
    atomsP_trace_length = length(md$atomsP %||% numeric()),
    equilibrationSnapshots = length(md$equilibrationSnapshotsA %||% list()),
    samplingSnapshots = length(md$samplingSnapshotsA %||% list()),
    pumpStatDim = paste(dim(md$pumpStat %||% matrix(numeric(), nrow = 0, ncol = 0)), collapse = "x"),
    meanPatternAssignmentDim = paste(dim(md$meanPatternAssignment %||% matrix(numeric(), nrow = 0, ncol = 0)), collapse = "x")
  )
  write_metrics(metrics_path, metrics)
}, error = function(e) {
  status <<- "error"
  err_txt <<- conditionMessage(e)
})

if (!identical(status, "ok")) {
  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write_metrics(
    metrics_path,
    list(
      status = status,
      language = "R",
      package = "CoGAPS",
      nPatterns = opt$k,
      nIterations = opt$n_iter,
      seed = opt$seed,
      sparseOptimization = isTRUE(opt$use_sparse_opt),
      nThreads = opt$cogaps_threads,
      asynchronousUpdates = isTRUE(opt$asynchronous_updates),
      openmpSupport = compiledWithOpenMPSupport(),
      outputFrequency = opt$output_frequency,
      checkpointInterval = opt$checkpoint_interval,
      checkpointOutFile = basename(checkpoint_path),
      nSnapshots = opt$n_snapshots,
      snapshotPhase = opt$snapshot_phase,
      takePumpSamples = isTRUE(opt$take_pump_samples),
      runtime_sec = runtime_sec,
      error = err_txt
    )
  )
  stop(err_txt)
}
