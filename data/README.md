# Case Study Data

This folder contains cached files used by the first-draft Bio-OCS CoGAPS case study.

## `processed/`

- `preprocessed_cells_hvg3000.h5ad`: processed cells x genes AnnData object.
- `cogaps_input_genesxcells_hvg3000_float64.h5ad`: genes x cells CoGAPS input.
- `preprocess_config_hvg3000.json`: preprocessing configuration used to create the cached object.

## `results/`

Legacy Python result location retained for compatibility:

- `cogaps_K7_seed2_iter2000.h5ad`: saved chosen Python CoGAPS model.
- `cogaps_K7_seed2_iter2000.metrics.json`: chosen-model metadata.

## `results_python/`

- `cogaps_K7_seed2_iter2000.h5ad`: saved chosen Python CoGAPS model in the language-specific results folder.
- `cogaps_K7_seed2_iter2000.metrics.json`: Python chosen-model metadata.
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

- `cogaps_K7_seed2_iter2000.rds`: saved chosen R CoGAPS model (local-only; not intended for GitHub).
- `cogaps_K7_seed2_iter2000.metrics.json`: R chosen-model metadata.
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

The default learner path loads a saved chosen CoGAPS result. Running one local model is optional, and the full SLURM sweep is instructor/advanced material. When the large full-result files are unavailable, the lighter exported pattern artifacts above provide a render-safe fallback for most of the interpretation workflow in both languages.
