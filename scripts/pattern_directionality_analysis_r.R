#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(zellkonverter)
  library(CoGAPS)
  library(SummarizedExperiment)
  library(Matrix)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

lookup_count <- function(tbl, key) {
  if (key %in% names(tbl)) {
    as.integer(tbl[[key]])
  } else {
    0L
  }
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "."
ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
PREPROCESSED_H5AD <- file.path(ROOT, "data", "processed", "preprocessed_cells_hvg3000.h5ad")
CHOSEN_RESULT_RDS <- file.path(ROOT, "data", "results_r", "cogaps_K7_seed2_iter2000.rds")
OUT_DIR <- file.path(ROOT, "data", "results_r")

TOP_N_GENES <- 15
MIN_CELLS_PER_PSEUDOBULK <- 20
CPM_PSEUDOCOUNT <- 0.5
STRONG_ABS_LOG2FC <- 0.5
MODERATE_ABS_LOG2FC <- 0.25
HIGH_CONSISTENCY <- 0.75

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

classify_direction <- function(mean_log2fc, sign_consistency, n_pairs) {
  if (is.na(mean_log2fc) || n_pairs == 0) {
    return(list(direction = "not_tested", evidence_strength = "insufficient_pairs"))
  }

  if (mean_log2fc > 0) {
    direction <- "up_in_stim"
  } else if (mean_log2fc < 0) {
    direction <- "down_in_stim"
  } else {
    direction <- "mixed_or_unclear"
  }

  abs_fc <- abs(mean_log2fc)
  if (abs_fc >= STRONG_ABS_LOG2FC && sign_consistency >= HIGH_CONSISTENCY) {
    strength <- "strong"
  } else if (abs_fc >= MODERATE_ABS_LOG2FC && sign_consistency >= 0.625) {
    strength <- "moderate"
  } else if (abs_fc < 0.1 || sign_consistency < 0.625) {
    strength <- "weak_or_mixed"
  } else {
    strength <- "suggestive"
  }

  if (identical(direction, "mixed_or_unclear")) {
    strength <- "weak_or_mixed"
  }

  list(direction = direction, evidence_strength = strength)
}

extract_top_pattern_genes <- function(A, top_n = 15) {
  do.call(
    rbind,
    lapply(pattern_columns(colnames(A)), function(pattern) {
      ordered <- sort(A[, pattern], decreasing = TRUE)
      ordered <- head(ordered, top_n)
      data.frame(
        pattern = pattern,
        gene = names(ordered),
        rank = seq_along(ordered),
        cogaps_weight = unname(ordered),
        row.names = NULL
      )
    })
  )
}

get_counts_assay <- function(sce) {
  assay_names <- assayNames(sce)
  if ("counts" %in% assay_names) {
    return(assay(sce, "counts"))
  }
  if ("X" %in% assay_names) {
    warning("counts assay not found; falling back to X. Directionality is best interpreted when raw counts are available.")
    return(assay(sce, "X"))
  }
  stop("Could not find a counts or X assay in the preprocessed H5AD.")
}

sum_gene_counts <- function(mat, gene_idx, cell_idx) {
  block <- mat[gene_idx, cell_idx, drop = FALSE]
  if (inherits(block, "sparseMatrix")) {
    return(as.numeric(Matrix::rowSums(block)))
  }
  rowSums(as.matrix(block))
}

group_library_size <- function(mat, cell_idx) {
  block <- mat[, cell_idx, drop = FALSE]
  if (inherits(block, "sparseMatrix")) {
    return(as.numeric(sum(block)))
  }
  sum(as.matrix(block))
}

make_pseudobulk <- function(sce, genes, group_cols, min_cells) {
  counts_mat <- get_counts_assay(sce)
  gene_idx <- match(genes, rownames(sce))
  if (anyNA(gene_idx)) {
    stop("Some genes were not found in the processed SingleCellExperiment.")
  }

  obs <- as.data.frame(colData(sce))
  obs$.cell_index <- seq_len(ncol(sce))

  keys <- do.call(paste, c(obs[group_cols], sep = "\r"))
  groups <- split(seq_len(nrow(obs)), keys)

  rows <- list()
  row_i <- 1L
  for (idx in groups) {
    group_df <- obs[idx, , drop = FALSE]
    if (nrow(group_df) < min_cells) {
      next
    }

    cell_idx <- group_df$.cell_index
    gene_counts <- sum_gene_counts(counts_mat, gene_idx, cell_idx)
    row <- as.list(group_df[1, group_cols, drop = FALSE])
    row$n_cells <- nrow(group_df)
    row$library_size_total_counts <- group_library_size(counts_mat, cell_idx)
    row$library_size_pattern_genes <- sum(gene_counts)
    for (j in seq_along(genes)) {
      row[[genes[[j]]]] <- gene_counts[[j]]
    }
    rows[[row_i]] <- row
    row_i <- row_i + 1L
  }

  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

add_cpm <- function(pseudobulk, genes) {
  out <- pseudobulk
  for (gene in genes) {
    out[[paste0(gene, "__cpm")]] <- ifelse(
      out$library_size_total_counts > 0,
      (out[[gene]] / out$library_size_total_counts) * 1e6,
      0
    )
  }
  out
}

paired_directionality <- function(pseudobulk_cpm, genes, pair_cols, context_cols) {
  if (!nrow(pseudobulk_cpm)) {
    return(list(pair_df = data.frame(), summary_df = data.frame()))
  }

  split_keys <- do.call(paste, c(pseudobulk_cpm[pair_cols], sep = "\r"))
  grouped <- split(seq_len(nrow(pseudobulk_cpm)), split_keys)

  pair_rows <- list()
  pair_i <- 1L

  for (idx in grouped) {
    group <- pseudobulk_cpm[idx, , drop = FALSE]
    if (!all(c("ctrl", "stim") %in% group$condition)) {
      next
    }
    ctrl <- group[group$condition == "ctrl", , drop = FALSE][1, ]
    stim <- group[group$condition == "stim", , drop = FALSE][1, ]

    pair_info <- as.list(ctrl[1, pair_cols, drop = FALSE])
    for (gene in genes) {
      ctrl_cpm <- as.numeric(ctrl[[paste0(gene, "__cpm")]])
      stim_cpm <- as.numeric(stim[[paste0(gene, "__cpm")]])
      log2fc <- log2(stim_cpm + CPM_PSEUDOCOUNT) - log2(ctrl_cpm + CPM_PSEUDOCOUNT)
      pair_rows[[pair_i]] <- c(
        context_cols,
        pair_info,
        list(
          gene = gene,
          ctrl_cpm = ctrl_cpm,
          stim_cpm = stim_cpm,
          log2fc_stim_vs_ctrl = log2fc,
          ctrl_n_cells = as.integer(ctrl$n_cells),
          stim_n_cells = as.integer(stim$n_cells)
        )
      )
      pair_i <- pair_i + 1L
    }
  }

  if (!length(pair_rows)) {
    return(list(pair_df = data.frame(), summary_df = data.frame()))
  }

  pair_df <- do.call(rbind.data.frame, c(pair_rows, stringsAsFactors = FALSE))
  numeric_cols <- c("ctrl_cpm", "stim_cpm", "log2fc_stim_vs_ctrl", "ctrl_n_cells", "stim_n_cells")
  for (col in numeric_cols) {
    pair_df[[col]] <- as.numeric(pair_df[[col]])
  }

  summary_group_cols <- unique(c("gene", setdiff(pair_cols, "replicate")))
  summary_keys <- do.call(paste, c(pair_df[summary_group_cols], sep = "\r"))
  summary_groups <- split(seq_len(nrow(pair_df)), summary_keys)

  summary_rows <- list()
  summary_i <- 1L

  for (idx in summary_groups) {
    gene_df <- pair_df[idx, , drop = FALSE]
    group_info <- as.list(gene_df[1, summary_group_cols, drop = FALSE])
    vals <- as.numeric(gene_df$log2fc_stim_vs_ctrl)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      next
    }
    mean_log2fc <- mean(vals)
    median_log2fc <- median(vals)
    n_pairs <- length(vals)
    n_up <- sum(vals > 0)
    n_down <- sum(vals < 0)
    n_zero <- sum(vals == 0)
    sign_consistency <- max(n_up, n_down, n_zero) / n_pairs
    direction <- classify_direction(mean_log2fc, sign_consistency, n_pairs)
    pvalue <- if (n_pairs >= 2 && !all(abs(vals) < .Machine$double.eps)) {
      tryCatch(wilcox.test(vals, mu = 0, alternative = "two.sided", exact = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }

    summary_rows[[summary_i]] <- c(
      context_cols,
      group_info,
      list(
        mean_log2fc_stim_vs_ctrl = mean_log2fc,
        median_log2fc_stim_vs_ctrl = median_log2fc,
        n_pairs = n_pairs,
        n_up_pairs = n_up,
        n_down_pairs = n_down,
        n_zero_pairs = n_zero,
        sign_consistency = sign_consistency,
        wilcoxon_pvalue = pvalue,
        direction = direction$direction,
        evidence_strength = direction$evidence_strength
      )
    )
    summary_i <- summary_i + 1L
  }

  summary_df <- do.call(rbind.data.frame, c(summary_rows, stringsAsFactors = FALSE))
  numeric_summary_cols <- c(
    "mean_log2fc_stim_vs_ctrl", "median_log2fc_stim_vs_ctrl", "n_pairs", "n_up_pairs",
    "n_down_pairs", "n_zero_pairs", "sign_consistency", "wilcoxon_pvalue"
  )
  for (col in numeric_summary_cols) {
    summary_df[[col]] <- as.numeric(summary_df[[col]])
  }
  summary_df$wilcoxon_fdr <- p.adjust(summary_df$wilcoxon_pvalue, method = "BH")

  list(pair_df = pair_df, summary_df = summary_df)
}

make_pattern_direction_summary <- function(global_summary) {
  do.call(
    rbind,
    lapply(split(global_summary, global_summary$pattern), function(df) {
      direction_counts <- table(df$direction)
      strength_counts <- table(df$evidence_strength)
      up_genes <- head(df$gene[df$direction == "up_in_stim"][order(df$mean_log2fc_stim_vs_ctrl[df$direction == "up_in_stim"], decreasing = TRUE)], 5)
      down_genes <- head(df$gene[df$direction == "down_in_stim"][order(df$mean_log2fc_stim_vs_ctrl[df$direction == "down_in_stim"], decreasing = FALSE)], 5)
      data.frame(
        pattern = df$pattern[[1]],
        n_top_genes = nrow(df),
        n_up_in_stim = lookup_count(direction_counts, "up_in_stim"),
        n_down_in_stim = lookup_count(direction_counts, "down_in_stim"),
        n_mixed_or_unclear = lookup_count(direction_counts, "mixed_or_unclear"),
        n_strong = lookup_count(strength_counts, "strong"),
        n_moderate = lookup_count(strength_counts, "moderate"),
        n_suggestive = lookup_count(strength_counts, "suggestive"),
        n_weak_or_mixed = lookup_count(strength_counts, "weak_or_mixed"),
        median_mean_log2fc = median(df$mean_log2fc_stim_vs_ctrl),
        top_up_genes = paste(up_genes, collapse = ", "),
        top_down_genes = paste(down_genes, collapse = ", "),
        row.names = NULL
      )
    })
  )
}

make_pattern_heatmap <- function(global_summary, output_path) {
  ordered <- global_summary[order(global_summary$pattern, global_summary$rank), ]
  ordered$label <- paste0(ordered$pattern, ":", ordered$gene)
  p <- ggplot(ordered, aes(x = "stim vs ctrl\nmean paired log2FC", y = reorder(label, rev(seq_along(label))), fill = mean_log2fc_stim_vs_ctrl)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(x = NULL, y = NULL, title = "Directionality of top CoGAPS pattern genes", fill = "mean paired log2FC") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 6))
  ggsave(output_path, plot = p, width = 5, height = max(8, nrow(ordered) * 0.18), dpi = 200)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sce <- zellkonverter::readH5AD(PREPROCESSED_H5AD, reader = "R")
result <- readRDS(CHOSEN_RESULT_RDS)
A <- as.matrix(getFeatureLoadings(result))
A <- ensure_pattern_matrix_orientation(A, nrow(sce), rownames(sce))
if (is.null(colnames(A))) {
  colnames(A) <- sprintf("Pattern%d", seq_len(ncol(A)))
} else {
  colnames(A) <- canonicalize_pattern_names(colnames(A))
}

top_genes <- extract_top_pattern_genes(A, TOP_N_GENES)
genes <- unique(top_genes$gene)

global_pb <- add_cpm(
  make_pseudobulk(
    sce = sce,
    genes = genes,
    group_cols = c("replicate", "condition"),
    min_cells = MIN_CELLS_PER_PSEUDOBULK
  ),
  genes
)
global_direction <- paired_directionality(
  global_pb,
  genes = genes,
  pair_cols = c("replicate"),
  context_cols = list(contrast_scope = "global")
)
global_summary <- merge(top_genes, global_direction$summary_df, by = "gene", all.x = TRUE)
global_summary <- global_summary[order(global_summary$pattern, global_summary$rank), ]

by_cell_pb <- add_cpm(
  make_pseudobulk(
    sce = sce,
    genes = genes,
    group_cols = c("replicate", "condition", "cell_type"),
    min_cells = MIN_CELLS_PER_PSEUDOBULK
  ),
  genes
)
cell_direction <- paired_directionality(
  by_cell_pb,
  genes = genes,
  pair_cols = c("replicate", "cell_type"),
  context_cols = list(contrast_scope = "by_cell_type")
)
cell_summary <- merge(top_genes, cell_direction$summary_df, by = "gene", all.x = TRUE)
cell_summary <- cell_summary[order(cell_summary$pattern, cell_summary$rank, cell_summary$cell_type), ]

pattern_summary <- make_pattern_direction_summary(global_summary)

write.csv(top_genes, file.path(OUT_DIR, "pattern_top_genes_directionality.csv"), row.names = FALSE)
write.csv(global_summary, file.path(OUT_DIR, "pattern_gene_directionality_global.csv"), row.names = FALSE)
write.csv(cell_summary, file.path(OUT_DIR, "pattern_gene_directionality_by_celltype.csv"), row.names = FALSE)
write.csv(pattern_summary, file.path(OUT_DIR, "pattern_direction_summary.csv"), row.names = FALSE)
write.csv(global_direction$pair_df, file.path(OUT_DIR, "pattern_gene_directionality_pairs_global.csv"), row.names = FALSE)
write.csv(cell_direction$pair_df, file.path(OUT_DIR, "pattern_gene_directionality_pairs_by_celltype.csv"), row.names = FALSE)

if (nrow(global_pb)) {
  write.csv(global_pb, file.path(OUT_DIR, "pseudobulk_counts_global.csv"), row.names = FALSE)
}
if (nrow(by_cell_pb)) {
  write.csv(by_cell_pb, file.path(OUT_DIR, "pseudobulk_counts_by_celltype.csv"), row.names = FALSE)
}

make_pattern_heatmap(global_summary, file.path(OUT_DIR, "pattern_gene_directionality_heatmap.png"))

cat("Wrote R directionality outputs to:", OUT_DIR, "\n")
