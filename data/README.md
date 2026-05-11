# Case Study Data

This folder contains cached files used by the first-draft Bio-OCS CoGAPS case study.

## `processed/`

- `preprocessed_cells_hvg3000.h5ad`: processed cells x genes AnnData object.
- `cogaps_input_genesxcells_hvg3000_float64.h5ad`: optional genes x cells CoGAPS input for full local Python/PyCoGAPS runs.
- `preprocess_config_hvg3000.json`: preprocessing configuration used to create the cached object.

## `results_r_k6_sparse_mt_t4_heavy/`

Selected R K6 sparse-on result location used by the default learner path:

- `cogaps_K6_seed2_iter2000.rds`: optional full R CoGAPS model object for local full-result workflows.
- `cogaps_K6_seed2_iter2000.metrics.json`: selected R model metadata.
- `pattern_gene_weights.tsv.gz`: lightweight exported gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: lightweight exported cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern from the selected K6 model.
- `pattern_correlations.csv`: correlations between pattern activity and the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact per-pattern cell-type summaries used in lightweight rendering.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries for optional extensions.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row summary per pattern.
- `pattern_top_genes_directionality.csv`: narrower top-gene subset used in directionality summaries.
- `pattern_gene_directionality_pairs_global.csv`: donor-paired global log2FC table.
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-paired cell-type log2FC table.
- `pattern_gene_directionality_heatmap.png`: render-ready heatmap summarizing top pattern-gene directionality.
- `pseudobulk_counts_global.csv`: global pseudobulk counts used to derive directionality summaries.
- `pseudobulk_counts_by_celltype.csv`: cell-type pseudobulk counts used to derive directionality summaries.

The full `.rds` file is optional for ordinary GitHub/classroom rendering. The lightweight files above are the cooking-show artifacts used by the main case-study text when the full model object is absent.

## `results_selected/`

Selected K6 Python/PyCoGAPS-compatible result location used by the secondary Python implementation:

- `cogaps_K6_seed2_iter2000.h5ad`: optional full Python saved model object for local full-result workflows.
- `cogaps_K6_seed2_iter2000.metrics.json`: selected Python model metadata.
- `pattern_gene_weights.tsv.gz`: lightweight exported Python gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: lightweight exported Python cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level Python pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern from the selected Python K6 model.
- `pattern_correlations.csv`: correlations between Python pattern activity and the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact Python per-pattern cell-type summaries used in lightweight rendering.
- `pattern_activity_by_replicate_condition.csv`: donor-aware Python activity summaries for optional extensions.
- `pattern_summary.csv`: compact Python per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for Python top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for Python top pattern genes.
- `pattern_direction_summary.csv`: one-row Python summary per pattern.
- `pattern_gene_directionality_pairs_global.csv`: donor-paired global log2FC table for the Python directionality analysis.
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-paired cell-type log2FC table for the Python directionality analysis.
- `pattern_gene_directionality_heatmap.png`: render-ready heatmap summarizing Python directionality estimates for top pattern genes.
- `pseudobulk_counts_global.csv`: global pseudobulk counts used to derive Python directionality summaries.
- `pseudobulk_counts_by_celltype.csv`: cell-type pseudobulk counts used to derive Python directionality summaries.

The full Python `.h5ad` file is optional for ordinary GitHub/classroom rendering. The lightweight files above are the secondary Python cooking-show artifacts used when the full Python model object is absent.

## `results/`

Legacy K7 Python result location retained for compatibility/reference:

- `cogaps_K7_seed2_iter2000.h5ad`: saved older K7 reference CoGAPS model.
- `cogaps_K7_seed2_iter2000.metrics.json`: older K7 reference metadata.

## `results_python/`

- Older K7/Python reference artifacts retained for compatibility and sensitivity checks.
- `cogaps_K7_seed2_iter2000.h5ad`: saved older K7 Python CoGAPS model in the language-specific results folder.
- `cogaps_K7_seed2_iter2000.metrics.json`: older K7 Python metadata.
- `summary_by_K.csv`: sweep summary used to justify the chosen rank.
- `per_run_metrics.csv`: run-level sweep metrics.
- `pattern_gene_weights.tsv.gz`: lightweight exported Python gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: lightweight exported Python cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level Python pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern from the chosen Python model.
- `pattern_correlations.csv`: correlations between Python pattern activity and the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact per-pattern cell-type summaries used in lightweight rendering.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries for optional extensions.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row summary per pattern.

## `results_r/`

- Older K7/R reference artifacts retained for compatibility and sensitivity checks.
- `cogaps_K7_seed2_iter2000.rds`: saved older K7 R CoGAPS model (local-only; not intended for GitHub).
- `cogaps_K7_seed2_iter2000.metrics.json`: older K7 R metadata.
- `pattern_gene_weights.tsv.gz`: lightweight exported R gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: lightweight exported R cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: R cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern from the chosen R model.
- `pattern_correlations.csv`: correlations between R pattern activity and the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact R per-pattern cell-type summaries used in lightweight rendering.
- `pattern_activity_by_replicate_condition.csv`: donor-aware R activity summaries for optional extensions.
- `pattern_summary.csv`: compact per-pattern interpretation table for the R implementation.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes from the R result.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes from the R result.
- `pattern_direction_summary.csv`: one-row summary per pattern from the R implementation.
- `pattern_top_genes_directionality.csv`: the narrower top-gene subset used in the R directionality analysis.
- `pattern_gene_directionality_pairs_global.csv`: donor-paired global log2FC table for the R directionality analysis.
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-paired cell-type log2FC table for the R directionality analysis.
- `pattern_gene_directionality_heatmap.png`: render-ready heatmap summarizing R directionality estimates for the top pattern genes.
- `pseudobulk_counts_global.csv`: global pseudobulk counts used to derive R directionality summaries.
- `pseudobulk_counts_by_celltype.csv`: cell-type pseudobulk counts used to derive R directionality summaries.

## `model_selection/`

Lightweight CSVs summarizing the lower-K, longer-iteration, high-K, and candidate-pattern model-selection analyses.

The default learner path uses the selected R K6 sparse-on cooking-show artifacts. Running one local model is optional, and the full SLURM sweep is instructor/advanced material. When full-result files are unavailable, the lighter exported pattern artifacts above provide a render-safe fallback for most of the interpretation workflow.
