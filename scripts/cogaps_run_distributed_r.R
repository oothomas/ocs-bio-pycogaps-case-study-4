#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
  library(BiocParallel)
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

choose_strata_cols <- function(obs, requested_cols) {
  requested_cols <- trimws(unlist(strsplit(requested_cols, ",", fixed = TRUE)))
  requested_cols <- requested_cols[nzchar(requested_cols)]
  available <- requested_cols[requested_cols %in% colnames(obs)]
  if (length(available) == 0) {
    return(character())
  }
  available
}

round_robin_assign <- function(index_groups, n_sets, seed) {
  sets <- vector("list", n_sets)
  for (i in seq_len(n_sets)) {
    sets[[i]] <- integer()
  }

  set.seed(seed)
  for (group_name in names(index_groups)) {
    group_idx <- index_groups[[group_name]]
    idx <- if (length(group_idx) > 1) {
      sample(group_idx, size = length(group_idx), replace = FALSE)
    } else {
      group_idx
    }
    target_sets <- rep(seq_len(n_sets), length.out = length(idx))
    for (i in seq_along(idx)) {
      sets[[target_sets[[i]]]] <- c(sets[[target_sets[[i]]]], idx[[i]])
    }
  }

  lapply(sets, sort)
}

build_stratified_cell_sets <- function(obs, n_sets, seed, strata_cols) {
  n_cells <- nrow(obs)
  if (length(strata_cols) == 0) {
    strata <- rep("all_cells", n_cells)
  } else {
    strata_parts <- lapply(strata_cols, function(col) as.character(obs[[col]]))
    strata_parts <- lapply(strata_parts, function(x) ifelse(is.na(x) | !nzchar(x), "NA", x))
    strata <- do.call(paste, c(strata_parts, sep = "||"))
  }

  indices_by_stratum <- split(seq_len(n_cells), strata)
  round_robin_assign(indices_by_stratum, n_sets = n_sets, seed = seed)
}

build_balanced_feature_sets <- function(mat, n_sets, seed) {
  feature_mean <- rowMeans(mat)
  feature_var <- apply(mat, 1, stats::var)
  ranking <- order(feature_var, feature_mean, decreasing = TRUE, na.last = TRUE)
  blocks <- split(ranking, ceiling(seq_along(ranking) / n_sets))
  indexed_blocks <- lapply(blocks, sort)
  round_robin_assign(indexed_blocks, n_sets = n_sets, seed = seed)
}

subset_summary_df <- function(obs, sets, strata_cols) {
  rows <- vector("list", length(sets))
  for (i in seq_along(sets)) {
    idx <- sets[[i]]
    sub_obs <- obs[idx, , drop = FALSE]
    base <- data.frame(
      set_id = sprintf("set%02d", i),
      n_cells = nrow(sub_obs),
      stringsAsFactors = FALSE
    )

    if (length(strata_cols) == 0) {
      base$stratum <- "all_cells"
      base$n_stratum <- nrow(sub_obs)
      rows[[i]] <- base
    } else {
      sub_obs$stratum <- do.call(
        paste,
        c(lapply(strata_cols, function(col) as.character(sub_obs[[col]])), sep = "||")
      )
      tab <- sort(table(sub_obs$stratum), decreasing = TRUE)
      rows[[i]] <- data.frame(
        set_id = sprintf("set%02d", i),
        n_cells = nrow(sub_obs),
        stratum = names(tab),
        n_stratum = as.integer(tab),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

feature_subset_summary_df <- function(mat, sets, feature_names) {
  rows <- lapply(seq_along(sets), function(i) {
    idx <- sets[[i]]
    feature_mean <- rowMeans(mat[idx, , drop = FALSE])
    feature_var <- apply(mat[idx, , drop = FALSE], 1, stats::var)
    data.frame(
      set_id = sprintf("set%02d", i),
      n_features = length(idx),
      mean_feature_mean = mean(feature_mean),
      median_feature_mean = stats::median(feature_mean),
      mean_feature_var = mean(feature_var),
      median_feature_var = stats::median(feature_var),
      top_feature = feature_names[idx[[which.max(feature_mean)]]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

extract_worker_summary <- function(first_pass_results, subset_sizes) {
  if (is.null(first_pass_results) || length(first_pass_results) == 0) {
    return(data.frame())
  }

  rows <- lapply(seq_along(first_pass_results), function(i) {
    res <- first_pass_results[[i]]
    md <- slot(res, "metadata")
    data.frame(
      worker = i,
      subset_size = subset_sizes[[i]],
      totalRunningTime = md$totalRunningTime %||% NA_real_,
      totalUpdates = md$totalUpdates %||% NA_real_,
      meanChiSq = md$meanChiSq %||% NA_real_,
      chisq_trace_length = length(md$chisq %||% numeric()),
      atomsA_trace_length = length(md$atomsA %||% numeric()),
      atomsP_trace_length = length(md$atomsP %||% numeric()),
      final_chisq = if (length(md$chisq %||% numeric()) > 0) tail(md$chisq, 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

extract_consensus_summary <- function(corr_to_mean_pattern) {
  if (is.null(corr_to_mean_pattern) || length(corr_to_mean_pattern) == 0) {
    return(data.frame())
  }

  rows <- lapply(seq_along(corr_to_mean_pattern), function(i) {
    vals <- corr_to_mean_pattern[[i]]
    data.frame(
      cluster = i,
      n_members = length(vals),
      min_corr_to_consensus = min(vals),
      mean_corr_to_consensus = mean(vals),
      max_corr_to_consensus = max(vals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

option_list <- list(
  make_option("--preprocessed-h5ad", type = "character", dest = "preprocessed_h5ad", help = "Processed cells x genes AnnData (.h5ad)"),
  make_option("--outdir", type = "character", default = "data/results_r_distributed", help = "Output directory [default %default]"),
  make_option("--k", type = "integer", help = "Number of patterns"),
  make_option("--seed", type = "integer", help = "Random seed"),
  make_option("--n-iter", type = "integer", dest = "n_iter", help = "Number of iterations"),
  make_option("--top-genes", type = "integer", dest = "top_genes", default = 50, help = "Top genes used for summary metrics [default %default]"),
  make_option("--stim-label", type = "character", dest = "stim_label", default = "stim", help = "Stimulated condition label [default %default]"),
  make_option("--distributed-mode", type = "character", dest = "distributed_mode", default = "single-cell", help = "Distributed mode: single-cell or genome-wide [default %default]"),
  make_option("--n-sets", type = "integer", dest = "n_sets", default = 4, help = "Number of distributed subsets [default %default]"),
  make_option("--cut", type = "integer", default = NA, help = "Pattern matching cut value [default K]"),
  make_option("--min-ns", type = "integer", dest = "min_ns", default = NA, help = "Minimum cluster size for pattern matching [default ceil(nSets/2)]"),
  make_option("--max-ns", type = "integer", dest = "max_ns", default = NA, help = "Maximum cluster size for pattern matching [default minNS + nSets]"),
  make_option("--bp-workers", type = "integer", dest = "bp_workers", default = NA, help = "BiocParallel workers [default nSets]"),
  make_option("--strata-cols", type = "character", dest = "strata_cols", default = "cell_type,condition,replicate", help = "Comma-separated metadata columns used to stratify cell subsets [default %default]"),
  make_option("--explicit-sets-rds", type = "character", dest = "explicit_sets_rds", default = NULL, help = "Optional RDS file containing explicit subsets to reuse [default create new sets]"),
  make_option("--output-frequency", type = "integer", dest = "output_frequency", default = 100, help = "Trace output frequency passed to CoGAPS [default %default]"),
  make_option("--quiet", action = "store_false", dest = "messages", default = TRUE, help = "Suppress CoGAPS status messages"),
  make_option("--write-sets-only", action = "store_true", dest = "write_sets_only", default = FALSE, help = "Write explicit sets and subset summaries, then exit before running CoGAPS"),
  make_option("--use-sparse-opt", action = "store_true", dest = "use_sparse_opt", default = TRUE, help = "Enable sparseOptimization [default]"),
  make_option("--no-sparse-opt", action = "store_false", dest = "use_sparse_opt", help = "Disable sparseOptimization"),
  make_option("--force-rerun", action = "store_true", default = FALSE, help = "Ignore successful cached outputs and rerun")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$preprocessed_h5ad) || is.null(opt$k) || is.null(opt$seed) || is.null(opt$n_iter)) {
  stop("Missing required arguments. Use --help for details.")
}

if (!(opt$distributed_mode %in% c("single-cell", "genome-wide"))) {
  stop("--distributed-mode must be either 'single-cell' or 'genome-wide'.")
}

outdir <- normalizePath(opt$outdir, winslash = "/", mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tag <- sprintf("K%d_seed%d_iter%d", opt$k, opt$seed, opt$n_iter)
run_label <- sprintf("%s_%s", tag, gsub("-", "_", opt$distributed_mode))
result_path <- file.path(outdir, sprintf("cogaps_%s.rds", run_label))
metrics_path <- file.path(outdir, sprintf("cogaps_%s.metrics.json", run_label))
diagnostics_rds_path <- file.path(outdir, sprintf("cogaps_%s.distributed_diagnostics.rds", run_label))
subset_summary_path <- file.path(outdir, sprintf("cogaps_%s.subset_summary.csv", run_label))
worker_summary_path <- file.path(outdir, sprintf("cogaps_%s.first_pass_worker_summary.csv", run_label))
consensus_summary_path <- file.path(outdir, sprintf("cogaps_%s.consensus_summary.csv", run_label))
explicit_sets_path <- file.path(outdir, sprintf("cogaps_%s.explicit_sets.rds", run_label))

if (file.exists(result_path) && file.exists(metrics_path) && !isTRUE(opt$force_rerun)) {
  prior_metrics <- tryCatch(fromJSON(metrics_path), error = function(e) NULL)
  if (!is.null(prior_metrics) && identical(prior_metrics$status, "ok")) {
    message(sprintf("[CACHE] status=ok; skipping %s", run_label))
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
  obs <- as.data.frame(colData(sce))
  feature_names <- rownames(sce)

  strata_cols <- character()
  explicit_sets_source <- "generated"
  if (!is.null(opt$explicit_sets_rds) && nzchar(opt$explicit_sets_rds)) {
    explicit_sets <- readRDS(opt$explicit_sets_rds)
    explicit_sets_source <- normalizePath(opt$explicit_sets_rds, winslash = "/", mustWork = TRUE)
  } else if (identical(opt$distributed_mode, "single-cell")) {
    strata_cols <- choose_strata_cols(obs, opt$strata_cols)
    explicit_sets <- build_stratified_cell_sets(
      obs = obs,
      n_sets = opt$n_sets,
      seed = opt$seed,
      strata_cols = strata_cols
    )
  } else {
    explicit_sets <- build_balanced_feature_sets(
      mat = data_matrix,
      n_sets = opt$n_sets,
      seed = opt$seed
    )
  }

  if (length(explicit_sets) != opt$n_sets) {
    stop("Number of explicit sets does not match --n-sets.")
  }
  saveRDS(explicit_sets, explicit_sets_path)

  if (identical(opt$distributed_mode, "single-cell")) {
    subset_summary <- subset_summary_df(obs, explicit_sets, strata_cols)
  } else {
    subset_summary <- feature_subset_summary_df(data_matrix, explicit_sets, feature_names = feature_names)
  }
  write.csv(subset_summary, subset_summary_path, row.names = FALSE)

  if (isTRUE(opt$write_sets_only)) {
    write_metrics(
      metrics_path,
      list(
        status = "sets_only",
        language = "R",
        package = "CoGAPS",
        distributedMode = opt$distributed_mode,
        nSets = opt$n_sets,
        seed = opt$seed,
        explicitSetsPath = basename(explicit_sets_path),
        explicitSetsSource = explicit_sets_source,
        strataCols = strata_cols,
        outputFrequency = opt$output_frequency
      )
    )
    quit(save = "no", status = 0)
  }

  params <- CogapsParams(
    nPatterns = opt$k,
    nIterations = opt$n_iter,
    seed = opt$seed,
    sparseOptimization = isTRUE(opt$use_sparse_opt)
  )
  params <- setParam(params, "distributed", opt$distributed_mode)

  cut_value <- if (is.na(opt$cut)) opt$k else opt$cut
  min_ns_value <- if (is.na(opt$min_ns)) ceiling(opt$n_sets / 2) else opt$min_ns
  max_ns_value <- if (is.na(opt$max_ns)) min_ns_value + opt$n_sets else opt$max_ns
  params <- setDistributedParams(
    params,
    nSets = opt$n_sets,
    cut = cut_value,
    minNS = min_ns_value,
    maxNS = max_ns_value
  )
  params <- setParam(params, "explicitSets", explicit_sets)

  bp_workers <- if (is.na(opt$bp_workers)) opt$n_sets else opt$bp_workers
  bp_param <- BiocParallel::MulticoreParam(workers = bp_workers, stop.on.error = TRUE)

  set.seed(opt$seed)
  result <- CoGAPS(
    data_matrix,
    params,
    nThreads = 1,
    messages = isTRUE(opt$messages),
    outputFrequency = opt$output_frequency,
    BPPARAM = bp_param
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
  if (!"condition" %in% colnames(obs)) {
    stop("The processed H5AD does not contain a 'condition' column in cell metadata.")
  }

  stim_indicator <- as.integer(as.character(obs$condition) == opt$stim_label)
  pattern_corrs <- vapply(
    pat_names,
    function(pattern) cor(sample_factors[, pattern], stim_indicator),
    numeric(1)
  )
  ifn_pattern <- names(which.max(pattern_corrs))
  ifn_corr <- unname(pattern_corrs[[ifn_pattern]])
  top_genes <- rownames(feature_loadings)[order(feature_loadings[, ifn_pattern], decreasing = TRUE)][seq_len(min(opt$top_genes, nrow(feature_loadings)))]

  md <- slot(result, "metadata")
  dx <- md
  saveRDS(dx, diagnostics_rds_path)

  worker_summary <- extract_worker_summary(dx$firstPass %||% list(), vapply(explicit_sets, length, integer(1)))
  if (nrow(worker_summary) > 0) {
    write.csv(worker_summary, worker_summary_path, row.names = FALSE)
  }

  consensus_summary <- extract_consensus_summary(dx$CorrToMeanPattern %||% list())
  if (nrow(consensus_summary) > 0) {
    write.csv(consensus_summary, consensus_summary_path, row.names = FALSE)
  }

  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  set_sizes <- vapply(explicit_sets, length, integer(1))

  metrics <- list(
    status = "ok",
    language = "R",
    package = "CoGAPS",
    result_path = basename(result_path),
    nPatterns = opt$k,
    nIterations = opt$n_iter,
    seed = opt$seed,
    sparseOptimization = isTRUE(opt$use_sparse_opt),
    distributedMode = opt$distributed_mode,
    nSets = opt$n_sets,
    cut = cut_value,
    minNS = min_ns_value,
    maxNS = max_ns_value,
    bpWorkers = bp_workers,
    explicitSetsPath = basename(explicit_sets_path),
    explicitSetsSource = explicit_sets_source,
    strataCols = strata_cols,
    subsetSizeMin = min(set_sizes),
    subsetSizeMedian = unname(stats::median(set_sizes)),
    subsetSizeMax = max(set_sizes),
    outputFrequency = opt$output_frequency,
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
    firstPassWorkers = length(dx$firstPass %||% list()),
    consensusClusters = nrow(consensus_summary),
    storedSubsets = length(dx$subsets %||% list())
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
      distributedMode = opt$distributed_mode,
      nSets = opt$n_sets,
      outputFrequency = opt$output_frequency,
      runtime_sec = runtime_sec,
      error = err_txt
    )
  )
  stop(err_txt)
}
