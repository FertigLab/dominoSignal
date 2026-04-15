# dominoSignal v1.6.0

## New Features

- Added `print()` and `show()` methods for `linkage_summary` objects to provide concise summary output.

## Bug Fixes

- Fixed `create_rl_map_cellphonedb()` handling of partner B complex mappings and gene assignment.
- Fixed `gene_network()` to avoid repeated prefixing of outgoing cluster names and to correctly subset outgoing signaling matrices.
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

- Added new `linkage_summary` class to summarize linkages in domino objects.
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

- Added chord plot of ligand expression targetting a specified receptor, with chord widths proportional to ligand expression per cell cluster.
- Added arguments to gene network plots to show communication between two clusters.
- Added filtering to signaling network plots to show outgoing signaling from specified clusters.

## Bug Fixes

- Fixed transcription factor-target linkages to exclude receptors within transcription factor regulon.
- Enabled `create_domino()` to run without providing a regulon list.
- Fixed ligand node sizing in gene network plots to correspond to the level of ligand expression.
