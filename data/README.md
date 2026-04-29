# Case Study Data

This folder contains cached files used by the first-draft Bio-OCS CoGAPS case study.

## `processed/`

- `preprocessed_cells_hvg3000.h5ad`: processed cells x genes AnnData object.
- `cogaps_input_genesxcells_hvg3000_float64.h5ad`: genes x cells CoGAPS input.
- `preprocess_config_hvg3000.json`: preprocessing configuration used to create the cached object.

## `results/`

- `cogaps_K7_seed2_iter2000.h5ad`: saved chosen CoGAPS model.
- `cogaps_K7_seed2_iter2000.metrics.json`: chosen-model metadata.
- `summary_by_K.csv`: sweep summary used to justify the chosen rank.
- `per_run_metrics.csv`: run-level sweep metrics.
- `pattern_gene_weights.tsv.gz`: lightweight exported gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: lightweight exported cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern from the chosen model.
- `pattern_correlations.csv`: correlations between pattern activity and the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact per-pattern cell-type summaries used in lightweight rendering.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries for optional extensions.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row summary per pattern.

The default learner path loads the saved chosen CoGAPS result. Running one local model is optional, and the full SLURM sweep is instructor/advanced material. When the large chosen-model file is unavailable, the lighter exported pattern artifacts above provide a render-safe fallback for most of the interpretation workflow.
