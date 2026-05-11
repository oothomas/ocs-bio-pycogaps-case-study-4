# Selected Python K6 CoGAPS Artifacts

This folder contains the secondary Python/PyCoGAPS-compatible artifacts for the selected `K = 6`, `seed = 2`, `n_iter = 2000`, sparse-on model.

The case-study text treats R/Bioconductor CoGAPS as the primary implementation. These Python artifacts are used as a secondary implementation and reproducibility/sensitivity check. Pattern labels are run-specific; the IFN-associated Python pattern is `Pattern5`, while the primary R K6 run labels the corresponding signal `Pattern4`.

## Full-Local File

- `cogaps_K6_seed2_iter2000.h5ad`: optional full Python saved model object for local full-result workflows.

This file is useful when a reviewer or developer wants to inspect the complete saved PyCoGAPS object or regenerate derived artifacts. It is intentionally optional for GitHub/classroom rendering because it is much larger than the lightweight cooking-show files.

## Lightweight Cooking-Show Files

These files are the default Python fallback when the full `.h5ad` is absent:

- `cogaps_K6_seed2_iter2000.metrics.json`: selected Python model metadata.
- `pattern_gene_weights.tsv.gz`: exported gene-weight (`A`) matrix.
- `pattern_cell_activities.tsv.gz`: exported cell-activity (`P`) matrix.
- `pattern_cell_activities_with_metadata.csv.gz`: cell-level pattern activities joined to metadata and embeddings.
- `pattern_top_genes.csv`: top genes per pattern.
- `pattern_correlations.csv`: pattern correlations with the stimulation indicator.
- `pattern_activity_by_celltype_condition.csv`: compact activity summaries by cell type and condition.
- `pattern_activity_by_replicate_condition.csv`: donor-aware activity summaries.
- `pattern_summary.csv`: compact per-pattern interpretation table.
- `pattern_gene_directionality_global.csv`: targeted global directionality for top pattern genes.
- `pattern_gene_directionality_by_celltype.csv`: targeted cell-type directionality for top pattern genes.
- `pattern_direction_summary.csv`: one-row directionality summary per pattern.
- `pattern_gene_directionality_pairs_global.csv`: donor-paired global log2FC values.
- `pattern_gene_directionality_pairs_by_celltype.csv`: donor-paired cell-type log2FC values.
- `pattern_gene_directionality_heatmap.png`: heatmap of global top-gene directionality.
- `pseudobulk_counts_global.csv`: global pseudobulk counts used for directionality summaries.
- `pseudobulk_counts_by_celltype.csv`: cell-type pseudobulk counts used for directionality summaries.

## Directionality Notes

CoGAPS weights identify genes associated with each pattern, but they do not by themselves indicate whether those genes are up- or down-regulated. Direction was estimated separately from the underlying count data by comparing IFN-beta stimulated cells to control cells with replicate-aware pseudobulk summaries.

- `Pattern5` is treated as the IFN-associated pattern for this secondary Python run.
- `Pattern5`: 14 of 15 top genes are up in stimulation by the global paired pseudobulk summary.
- This supports using directionality as evidence for interpreting the IFN-associated pattern, while keeping non-IFN or secondary pattern labels cautious.

## Caveats

- This is a targeted directionality check for pattern genes, not a genome-wide differential-expression analysis.
- The by-cell-type summaries are descriptive and depend on cell counts per donor/condition/cell-type stratum.
- The current dataset is an IFN-beta stimulation PBMC dataset; any dengue interpretation should be framed as biological relevance to IFN/inflammatory programs, not direct dengue infection evidence.
