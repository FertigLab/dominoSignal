# dominoSignal v1.7.1 (October 12, 2026)

## BREAKING CHANGES

Changes are noted as breaking when a call on the same inputs gives a different result or a call that used to work now errors. A number of bug fixes are included here (prefixed with [BF]), as the fixes alter the behavior of functions.

### Renamed or removed functions, arguments, and data

- Removed the `use_clusters` argument from `create_domino()`, as `TRUE` is the only option. The `counts`, `z_scores`, and `clusters` arguments are now required and no longer default to `NULL`.
- Renamed the `min_tf_pval` argument of `build_domino()` to `max_tf_pval`, as it as a maximum p-value threshold. Objects built with previous versions (which store `min_tf_pval` in `build_vars`) still pass validity checks.
- Renamed the first argument of `count_linkage()` and `test_differential_linkages()` from `linkage_summary` to `link_summary` to avoid conflict with `linkage_summary()` class.
- Removed the `test_name` argument from `test_differential_linkages()`, as Fisher's exact test is the only option.
- Renamed the `return` argument of `dom_network_items()` to `which_return` to avoid conflict with `return()` function.
- `obtain_circos_expression()` and `render_circos_ligand_receptor()` are no longer exported. Use `circos_ligand_receptor()`.
- Renamed elements of the `PBMC` example dataset: `RNA_count_tiny` is now `count_tiny` and `RNA_zscore_tiny` is now `zscore_tiny`.
- Removed `mock_linkage_summary()`. An example `linkage_summary()` object is now provided in the `LinkageSummary` dataset.

### Changed defaults

- Changed the default of `remove_rec_dropout` in `create_domino()` to `FALSE` to reflect current recommended usage.
- Changed the default of `remove_rec_dropout` in `cor_scatter()` to `NULL`, which uses the value stored in the `domino()` object by `create_domino()`. Objects created before v1.7.1 do not store this value and now default to `FALSE` (previously `TRUE`).
- The `plot_differential_linkages()` colors now fun from the `gradient` minimum color (red, now customizable) to the maximum color (grey, now customizable) regardless of `stat_ranking`. Previously, `stat_ranking = "descending"` colored the maximum value red.

### Stricter input checks and new errors

- `domino()` objects now have validity checks: cell-level slots must share the same cells in the same order, and stored parameters, `clust_de`, and `cor` must be consistent within the object. Objects created with earlier versions may fail `validObject()`, in which case they should be recreated.
- [BF] `linkage_summary()` objects now have validity checks: the first column of `subject_meta` and the names of `subject_linkages` must match `subject_names` in the same order to prevent mismatches when matching by position only. Objects created with earlier versions may fail `validObject()` if the order differed, in which case they should be recreated and results may change.
- Added argument validation to exported functions. Inputs that were previously accepted or failed later with unclear errors now error immediately with more informative errors.
- [BF] `create_domino()` matches cells across `counts`, `z_scores`, `features`, and `clusters` by name instead of only by position. `clusters` must be named by cell. Cells not shared by all inputs are dropped with a warning, and the function errors if no cells are shared. Results will differ for inputs whose cells were not in the same order.
- [BF] `create_domino()` excludes `rl_map` rows that do not pair one receptor with one ligand with a warning, and errors if no rows remain. Previously, `gene_B` was silently treated as the receptor if `gene_A` was a ligand.
- [BF] `summarize_linkages()` requires `domino()` objects built with `build_domino()` and subject names present in the first column of `subject_meta`. `subject_meta` is reordered to match `subject_names` instead of being matched by position. Results may change if order differed previously.
- `signaling_network()` now errors instead of warning and returning `NULL` when no signaling is found.
- `gene_network()` no errors instead of warning when the input `domino()` object has not been built.

### Changed results

- When `create_domino()` TF selection method is `all` or `variable`, cluster labels are kept and per-cluster networks are built using TFs ranked by maximum correlation with with receptors expressed in each cluster. Previously, clusters were set to an empty factor and a single `clust` list was returned, which was incompatible with downstream exploration and plotting functions.
- [BF] Receptor-TF correlations that are not calculated in `create_domino()` (receptors in the TF's regulon, or TFs with all-zero scores in the retained cells) are now `NA` instead of zero. Receptor complexes with such a component now have an `NA` median correlation and are no longer linked to that TF.
- [BF] `create_domino()` no longer drops receptors with non-syntactic names (such as those containing spaces or hyphens) from the signaling network.
- [BF] `build_domino()` keeps a receptor when it is the only receptor gene passing `min_rec_percentage` in a cluster. Previously, it was silently excluded from the network.
- [BF] `build_domino()` calculates signaling scores for clusters with only one incoming ligand. Previously, these scores were always zero.
- [BF] `build_domino()` averages the component genes of a complex ligand when it is the only valid ligand for a cluster, instead of summing the components as separate ligands. Complexes with missing components are dropped consistently.
- [BF] `create_rl_map_cellphonedb()` skips interactions whose components lack an ortholog during gene conversion, instead of keeping unconverted gene names. Names of single-protein partner B entries now replace spaces with underscores, matching partner A.
- [BF] `count_linkage()` and `test_differential_linkages()` now filter by `subject_names`. Previously, the argument was ignored.
- [BF] `count_linkage()` now counts the correct subjects when the subject name column of `subject_meta` is a factor.
- [BF] `summarize_linkages()` no longer records spurious `NA <- NA` linkages for subjects or clusters without linkages.
- `signaling_heatmap()` and `incoming_signaling_heatmap()` with `scale = "log"` now plot `log10(x + 1)` instead of `log10(x)` to avoid `log10(0)`.
- [BF] `gene_network()` no longer counts ligand expression more than once for vertex sizes when `OutgoingSignalingClust` is used with multiple receptor clusters.

## New Features

- `create_domino()` stores its parameters in the `domino()` object, and `create_domino()` and `build_domino()` record the dominoSignal version used. These are returned by `dom_info()` and reported by `print()` and `show()`, which now give the same output for `domino()` objects.
- `create_domino()` accepts a dense `matrix` or `data.frame` for `counts` and converts it to a sparse `dgCMatrix`.
- `build_domino()` messages when no TFs in a cluster pass the `max_tf_pval` threshold, explaining why no signaling is inferred.
- Added a [`subset`](../reference/subset-linkage_summary-method.html) method for `linkage_summary()` objects to filter by subject names or metadata.
- Added `gradient` argument to `plot_differential_linkages()` to set the colors of the test statistic.
- Added `dom_to_df()` function to create data.frame of signaling results from `domino()` object.
- Added example datasets `DominoObjects` (created and built `domino()` objects) and `LinkageSummary` (a `linkage_summary()` object and differential linkage results). Added `regulon_list_tiny` to `SCENIC` and `rl_map_tiny` to `CellPhoneDB`.

## Bug Fixes

- Fixed `create_domino()` failing with a custom `rl_map` that lacks `name_A` and `name_B` columns.
- Fixed `create_domino()` failing when a single feature is retained with `tf_selection_method = "variable"`, and added validation for `tf_variance_quantile`.
- Fixed `create_domino()` failing when receptor genes are missing from `z_scores`. These receptors are now excluded from correlation calculations with a warning.
- Error messages now refer to `create_domino()` and `build_domino()` instead of `domino_create` and `domino_build`.
- Fixed `create_rl_map_cellphonedb()` failing when receptor annotations are missing from the protein table.
- Fixed `create_rl_map_cellphonedb()` failing when no interactions remain after ortholog conversion. An empty `rl_map` is now returned with a warning.
- Fixed `test_differential_linkages()` failing when comparing more than two groups (the odds ratio is `NA` in this case).
- Fixed `feat_heatmap()` failing when `ann_cols = FALSE` and a title is used.
- Fixed `signaling_network()` producing `NaN` edge weights when normalizing clusters without signaling (0/0).
- `gene_network()` gives an informative error when `clust` is not provided.
- Fixed `gene_network()` failing when there are no linkages to plot.
- Fixed `plot_differential_linkages()` failing when the test statistic contains `NA` values (such as `odds.ratio`). `stat_range` is only restricted to [0, 1] for p-values.

## Documentation

- Updated vignettes for consistent terminology, corrected pySCENIC commands, and improved accessibility.
- Added vignette for the differential signaling workflow
- Regenerated example data with dominoSignal v1.7.1. Examples now load the prebuilt `DominoObjects` and `LinkageSummary` datasets instead of re-running `example(build_domino)`.
- Documented sources and contents of example datasets.
- Corrected documentation of `tf_variance_quantile` in `create_domino()` to state that the coefficient of variation is used and higher values keep fewer features.
- Fixed `create_rl_map_cellphonedb()` examples to use genes and a toy ortholog table that return results.


# dominoSignal v1.6.0 (May 5, 2026)

## New Features

- Added [`print`](../reference/print-linkage_summary-method.html) and [`show`](../reference/show-linkage_summary-method.html) methods for `linkage_summary()` objects to provide concise summary output.

## Bug Fixes

- Fixed `create_rl_map_cellphonedb()` handling of partner B complex mappings and gene assignment.
- Fixed `gene_network()` to avoid repeated prefixing of outgoing cluster names and to correctly subset outgoing signaling matrices.
- Fixed `gene_network()` to only include receptor and TF nodes that are associated with a ligand if `OutgoingSignalingClust` is used.
- Fixed `gene_network()` to accumulate ligand expression across clusters for ligand node scaling.
- Fixed `signaling_network()` to assign undefined (`NA`) vertex sizes to 0 when scaling by signaling.
- Fixed `dom_linkages()` with `by_cluster = TRUE` and `link_type = "tf-receptor"` to return `clust_tf_rec`.
- Fixed `dom_signaling(cluster = ...)` to return the selected cluster matrix via list indexing.

## Documentation

- Added figure alt text to images in vignettes for accessibility.
- Updated pkgdown and vignette links to use working URLs.
- Updated README/index documentation links and citation text to current release metadata.

# dominoSignal v1.4.1 (March 25, 2026)

## Other changes

- Updated maintainer information.

# dominoSignal v1.2.0 (April 15th, 2025)

## Bug Fixes

- Fixed `circos_ligand_receptor()` to not fail when rl_map includes ligands not present in the expression matrix. Missing ligands are excluded with informative message.
- Fixed `create_domino()` to prevent overwriting signaling matrix with `NULL` when `complexes = TRUE` but no complexes are found to have active signaling.

# dominoSignal v1.0.0 (November 6, 2024)

- Accepted to Bioconductor in release 3.20.

## Bug Fixes

- Disabled exact p-value computation for correlation test between receptor expression and features to prevent repeated warning messages due to inevitable tied ranks during Spearman correlation calculation in `create_domino()`.

## Documentation

- Updated vignette download instructions to use the Bioconductor URL.
- All vignettes explicitly state seed used when executing code if applicable.
- Example code runs with `echo = FALSE` to reduce output verbosity in documentation.
- `create_domino()` examples run with `verbose = FALSE` to reduce extensive output in documentation.
- Vignette regarding dominoSignal object structure explains the purpose of downloading and importing data with `BiocFileCache` to demonstrate applications on large real data objects.
- Fixed example code for `circos_ligand_receptor()` color customization and `cor_heatmap()` boolean representation.
- Updated non-functional links to correct URLs.

# dominoSignal v0.99.2-alpha (May 15, 2024)

## BREAKING CHANGES

- Renamed package from "domino2" to "dominoSignal".

## Documentation

- Updated vignettes to demonstrate pipeline on data formatted as `SingleCellExperiment` objects.
- Added SCENIC tutorial vignette in place of deprecated example scripts.

# dominoSignal v0.2.2-alpha (December 20, 2023)

## New Features

- Added new `linkage_summary()` class to summarize linkages in `domino()` objects.
- Added helper functions to count linkages and compare between `domino()` objects.
- Added plotting function for differential linkages.

# dominoSignal v0.2.1-alpha (November 3, 2023)

## New Features

### Function Inputs

- Standardized input formats for receptor-ligand databases, transcription factor activity scores, and regulon gene lists to support alternative databases and transcription factor activity inference methods.
- Added helper functions to reformat pySCENIC outputs and CellPhoneDB database files to standardized input formats.
- Added `host` option for gene ortholog conversions using `biomaRt` for access to maintained mirrors.

### Improved Linkages

- Implemented assessment of transcription factor linkages with heteromeric receptor complexes based on correlation between transcription factor activity and all receptor component genes
- Implemented assessment of complex ligand expression as the mean of component gene expression for plotting functions.
- Added minimum threshold parameter for the percentage of cells in a cluster expressing a receptor gene.
- Added linkage slots for active receptors per cluster, transcription factor-receptor linkages per cluster, and incoming ligands for active receptors within each cluster.

### Plotting Functions

- Added chord plot of ligand expression targeting a specified receptor, with chord widths proportional to ligand expression per cell cluster.
- Added arguments to gene network plots to show communication between two clusters.
- Added filtering to signaling network plots to show outgoing signaling from specified clusters.

## Bug Fixes

- Fixed transcription factor-target linkages to exclude receptors within transcription factor regulon.
- Enabled `create_domino()` to run without providing a regulon list.
- Fixed ligand node sizing in gene network plots to correspond to the level of ligand expression.
