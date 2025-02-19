# dominoSignal v1.0.5

- circos_ligand_receptor will not fail in cases where the rl_map of a domino objects includes ligands for a receptor where the ligands are not present in the expression matrix. A message is returned if ligands from the rl_map had to be excluded from the plot

# dominoSignal v1.0.4

- refactorization of circos_ligand_receptor function to remove use of grepl-based regular expressions to reformat cell and molecule names from a domino object into a data frame for plotting a receptor circos plot. Component functions for creating the data frame of ligand expression centered on a single receptor and to render this data frame as a circos plot are now included in the package as "obtain_circos_expression" and "render_circos_ligand_receptor", respectively.

# dominoSignal v1.0.3

- Added functions for calculation of mean gene expression among components of a complex using purrr functions


# dominoSignal v1.0.2

## Pkgdown Site Customization Scripts

- Restored _pkgdown.yml and index scripts that specify site building parameters for https://FertigLab.github.io/dominoSignal to improve formatting.

# dominoSignal v1.0.1

## Receptor Complex Bugfix

- Resolved issue where if create_domino was run with complexes=TRUE and no complexes were found to have active signaling, the full signaling matrix would be replaced by a NULL value.

## GitHub Actions

- Restoration of .github/workflows scripts for automatic build checks.

# dominoSignal v1.0.0

- Version number update to signify acceptance to bioconductor in release 3.20

# dominoSignal v0.99.4

## Vignettes

- Updated download instructions to use the bioconductor URL

# dominoSignal v0.99.3

## create_domino Function Warnings

- Disabled exact p-value computation for the correlation test between receptor expression and features to prevent repeated warning messages due to tied ranks during Spearman correlation calculation 

## Vignettes

- Updated non-functional links to functional URLs
- All vignettes explicitly state the seed used when executing their code
- The dominoSignal Object vignette states the purpose of the code used to download and import data from a BioFileCache to demonstrate applications of the domino objects on a real data object that is too large to include inside the package

## Testing Data

- Data used for examples and unit tests is now stored in the package's data/ directory.
- All usage of the triple colon operator to access internal data from sysdata have been replaced with usage of data()
- Linkage summaries for demonstrating differential linkage testing and plotting are generated within unit test scripts rather than being stored within the package.
- Examples are run with echo = FALSE to cut down on lines of code printed in example pages.
- Examples for create_domino are run with verbose = FALSE
- Fixed example of changing colors for circos_ligand_receptor
- Fixed example of boolean representation of cor_heatmap

# dominoSignal v0.99.2

## Package Name

- Update of package name from "domino2" to "dominoSignal"

## Bioconductor Standards

- Update to vignettes presenting application of the DominoSignal pipeline on data formatted as a SingleCellExperiment object
- Implemented caching of example data by BiocCache to meet package size limits
- Removal of deprecated scripts for running SCENIC. Tutorials for running SCENIC are still present in vignettes
- Corrected BiocCheck notes pertaining to coding practices including paste in conditional statements, functions with dontrun examples, usage of seq_len or seq_along in place of seq, and usage of vapply in place of sapply


# dominoSignal v0.99.1

- Update to Bioconductor version numbering conventions for package submission

# dominoSignal v0.2.2

## Linkage functions
- Addition of new class to summarize linkages in objects
- Addition of helper functions to count linkages and compare between objects
- Plotting function for differential linkages

## Package structure
- Adjustments made to meet Bioconductor standards

# dominoSignal v0.2.1

## Updates to domino object construction
- Uniform formats for inputs of receptor - ligand interaction databases, transcription factor activity features, and regulon gene lists for operability with alternative databases and transcription factor activation inference methods
- Helper functions for reformatting pySCENIC outputs and CellPhoneDB database files to domino-readable uniform formats
- Assessment of transcription factor linkage with receptors that function as a heteromeric complex based on correlation between transcription factor activity and all receptor component genes
- Assessment of complex ligand expression as the mean of component gene expression for plotting functions
- Minimum threshold for the percentage of cells in a cluster expressing a receptor gene for the receptor to be called active within the cluster
- Additional linkage slots for active receptors in each cluster, transcription factor - receptor linkages for each cluster, and incoming ligands for active receptors on each cluster

## Plotting functions
- Chord plot of ligand expression targeting a specified receptor where chord widths correspond to the quantity of ligand expression by each cell cluster
- Signaling networks showing only outgoing signaling from specified cell clusters
- Gene networks between two cell clusters

## Bugfixes
- Added host option for gene ortholog conversions using `{biomaRt}` for access to maintained mirrors
- Transcription factor - target linkages are now properly stored so that receptors in a transcription factor's regulon are excluded from linkage
- Ligand nodes sizes in gene networks correspond to quantity of ligand expression
- `create_domino()` can be run without providing a regulon list
- References to the host GitHub repository have been updated to [Elisseeff-Lab](https://github.com/Elisseeff-Lab/domino)
