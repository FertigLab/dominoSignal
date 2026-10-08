# dominoSignal v1.7.1 (in development)

## BREAKING CHANGES

- Updated name of first argument of `count_linkage()` and `test_differential_linkages()` functions to `link_summary` to avoid conflict with `linkage_summary()` class.
- Removed `test_name` argument from `test_differential_linkages()` function, as Fisher's Exact Test is the only option.
- Default for `remove_rec_dropout` in `create_domino()` has been changed to `FALSE` to reflect current recommended usage.
- Argument `use_clusters` removed from `create_domino()`, as `TRUE` is the only option.
- `domino()` validity check now requires that cell inputs have aligned names. Previously built objects (before the bug fix for cluster matching by positions) will now fail the validity check and should be re-run.
- When `create_domino()` TF selection method is set to `all` or `variable`, per cluster networks (using TFs ranked by maximum correlation in expressed receptors) are generated rather than the previous `clust` list. This ensures compatibility with downstream exploration and visualization functions but will change returned results compared to previous versions.
- Interactions in `rl_map` that do not pair one receptor (`R`) with one ligand (`L`) are excluded in `create_domino()` with a warning and the function errors if no rows remain (previously gene_B was silently treated as a receptor in those rows).
- Component genes of a complex ligand when it is the only valid ligand for a cluster are now averaged in `build_domino()` instead of summing component genes as separate rows. Complexes with missing components drop consistently. Signaling scores for affected clusters will differ from previous versions.
- The `summarize_linkages()` function requires `domino()` objects built with `build_domino()` and subject names present in the first column of `subject_meta`. The function reorders `subject_meta` to match `subject_names` to avoid issues with positional matching.
- Validity check for `linkage_summary()` object now requires the first column of `subject_meta` and the names of `subject_linkages` to match `subject_names` in the same order.

## New Features

- Added `dom_to_df()` function to create data.frame of signaling results from domino object.
- Added [`subset`](../reference/subset-linkage_summary-method.html) for `linkage_summary()` class to allow for filtering of objects by subject names or metadata.
- Added British-spelling synonyms for relevant functions (`summarise_linkages()`, `dom_signalling()`, `signalling_heatmap()`, `incoming_signalling_heatmap()`, `signalling_network()`).
- Added `gradient` argument to `plot_differential_linkages()` function for statistic coloring.
- Parameters for `domino()` object creation (used in `create_domino()`) are now stored in object
- Dense `matrix` or `data.frame` accepted for `counts` input to `create_domino()` and converted to sparse `dgCMatrix`.
- Default for `remove_rec_dropout` in `cor_scatter()` has been changed to `NULL` and will check `domino()` metadata for creation of the object (otherwise will be set to `FALSE`).

## Bug Fixes

- Fixed `gene_network()` function skip iterating on linkages if length is zero.
- Fixed `create_domino()` to work with custom rl_map without names_A and names_B columns
- Fixed `build_domino()` to retain receptor gene names when only one receptor passes expression threshold in a cluster, rather than silently excluding it from signaling network.
- Fixed `build_domino()` to compute signaling scores for clusters with only one incoming ligand instead of always returning zero.
- Fixed `count_linkage()` and `test_differential_linkages()` functions to use `subject_names` argument to filter `linkage_summary()` object before calculating results.
- Fixed `summarize_linkages()` iteration over linkage pairs (receptor-TF or ligand-receptor) to skip empty cases (removing spurious NA <- NA results when no links are present).
- Fixed `plot_differential_linkages()` to accept statistic ranges outside of [0, 1] for odds.ratio.
- Fixed `create_domino()` setting clusters to empty factor if using `tf_selection_method` of `variable` or `all`.
- Removed default value of `NULL` for required arguments of `counts`, `zscores`, and `clusters` in `create_domino()`.
- Fixed cluster matching to cells by position only by aligning cell inputs by name in `create_domino()`.
- The maximum p-value threshold for cluster based selection of TFs in `build_domino()` has argument name `max_tf_pval` instead of incorrect `min_tf_pval`. (The previous `min_tf_pval` name stored in `build_vars` of the `domino()` object should still pass validity checks.)
- Fixed `create_domino()` dropping receptors with non-syntactic names (such as those containing spaces or hyphens) from signaling network which were renamed in correlation matrix.
- Fixed `create_domino()` failing when single feature is retained by `tf_selection_method` is set to `variable` and added validation for `tf_variance_quantile` argument.
- Fixed `count_linkage()` counting the wrong subjects when subject name column of `subject_meta` is a factor.
- Fixed `create_rl_map_cellphonedb()` keeping components that lack an ortholog (resulting in unconverted gene names) instead of skipping interactions
- Fixed `create_rl_map_cellphonedb()` to handle missing receptor annotations in protein table input instead of crashing.
- When `signaling_heatmap()` and `incoming_signaling_heatmap()` are called with `scale = "log"` the transformation is now `log10(x+1)` instead of `log10()` to avoid `log10(0)` issues.
- `clust` is now a required argument of `gene_network()` (previously defaulted to `NULL`, which produced no plot).
- Fixed `gene_network()` counting ligand expression more than once for vertex sizes when `OutgoingSignalingClust` is used with multiple receptor clusters.

## Documentation

- Added vignette for differential signaling workflow.
- Example data regenerated with `dominoSignal` version 1.7.1.
- Corrected documentation of `tf_variance_quantile` in `create_domino()` to indicate coefficient of variation is used and higher values keep fewer features.
- Error message referring to `domino_create` instead of `create_domino()` in `build_domino()` has been corrected.
- When no TFs pass `max_tf_pval` threshold in `build_domino()`, users now receive a message explaining why no signaling is inferred.
- Fixed `create_rl_map_cellphonedb()` examples to use genes and toy ortholog table that returns results.

# dominoSignal v1.6.0

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

# dominoSignal v1.4.1

- Updated maintainer information.

# dominoSignal v1.2.0

## Bug Fixes

- Fixed `circos_ligand_receptor()` to not fail when rl_map includes ligands not present in the expression matrix. Missing ligands are excluded with informative message.
- Fixed `create_domino()` to prevent overwriting signaling matrix with `NULL` when `complexes = TRUE` but no complexes are found to have active signaling.

# dominoSignal v1.0.0

- Accepted to Bioconductor in release 3.20.

## Bug Fixes

- Disabled exact p-value computation for correlation test between receptor expression and features to prevent repeated warning messages due to inevitable tied ranks during Spearman correlation calculation in `create_domino()`.

## Documentation

- Updated vignette download instructions to use the Bioconductor URL
- All vignettes explicitly state seed used when executing code if applicable.
- Example code runs with `echo = FALSE` to reduce output verbosity in documentation
- `create_domino()` examples run with `verbose = FALSE` to reduce extensive output in documentation.
- Vignette regarding dominoSignal object structure explains the purpose of downloading and importing data with `BiocFileCache` to demonstrate applications on large real data objects.
- Fixed example code for `circos_ligand_receptor()` color customization and `cor_heatmap()` boolean representation.
- Updated non-functional links to correct URLs.

# dominoSignal v0.99.2-alpha

- Package renamed from "domino2" to "dominoSignal".

## Documentation

- Updated vignettes to demonstrate pipeline on data formatted as `SingleCellExperiment` objects.
- Added SCENIC tutorial vignette in place of deprecated example scripts

# dominoSignal v0.2.2-alpha

## New Features

- Added new `linkage_summary()` class to summarize linkages in domino objects.
- Added helper functions to count linkages and compare between domino objects.
- Added plotting function for differential linkages.

# dominoSignal v0.2.1-alpha

## New Features

### Function Inputs

- Standardized input formats for receptor-ligand databases, transcription factor activity scores, and regulon gene lists to support alternative databases and transcription factor activation inference methods.
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
