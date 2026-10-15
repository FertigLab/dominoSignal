# Create a [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md) object and prepare it for network construction

This function reads in a receptor ligand signaling database, cell level
features of some kind (ie. output from pySCENIC), z-scored single cell
data, and cluster id for single cell data, calculates a correlation
matrix between receptors and other features (this is transcription
factor module scores if using pySCENIC), and finds features enriched by
cluster. It will return a
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object prepared for
[`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md),
which will calculate a signaling network.

## Usage

``` r
create_domino(
  rl_map,
  features,
  counts,
  z_scores,
  clusters,
  tf_targets = NULL,
  verbose = TRUE,
  use_complexes = TRUE,
  rec_min_thresh = 0.025,
  remove_rec_dropout = FALSE,
  tf_selection_method = "clusters",
  tf_variance_quantile = 0.5
)
```

## Arguments

- rl_map:

  Data frame where each row describes a receptor-ligand interaction with
  required columns gene_A & gene_B including the gene names for the
  receptor and ligand and type_A & type_B annotating if genes A and B
  are a ligand (L) or receptor (R). Each interaction must pair one
  receptor with one ligand (rows that do not, such as ligand-ligand or
  missing types) are excluded with a warning.

- features:

  Either a path to a csv containing cell level features of interest (ie.
  the auc matrix from pySCENIC) or named matrix with cells as columns
  and features as rows.

- counts:

  Counts matrix for the data. Used to threshold receptors based on
  expression (`rec_min_thresh`), to calculate percentage of clusters
  with receptor expression, and to remove cells with zero receptor
  counts if `remove_rec_dropout` is TRUE. Dense matrices are converted
  to sparse `dgCMatrix` for storage.

- z_scores:

  Matrix containing z-scored expression data for all cells with cells as
  columns and features as rows. Receptor genes that pass
  `rec_min_thresh` but are not present in the z_scores matrix are
  excluded from correlation calculations with a warning.

- clusters:

  Named factor containing cell cluster with names as cells.

- tf_targets:

  Optional. A list where names are transcription factors and the stored
  values are character vectors of genes in the transcription factor's
  regulon.

- verbose:

  Boolean indicating whether or not to print progress during
  computation.

- use_complexes:

  Boolean indicating whether you wish to use receptor/ligand complexes
  in the receptor ligand signaling database. If FALSE, receptor/ligand
  pairs where either functions as a protein complex will not be
  considered when constructing the signaling network.

- rec_min_thresh:

  Minimum expression level of receptors by cell. Default is 0.025 or 2.5
  percent of all cells in the data set. This is important when
  calculating correlation to connect receptors to transcription
  activity. If this threshold is too low then correlation calculations
  will proceed with very few cells with non-zero expression.

- remove_rec_dropout:

  Whether to remove receptors with 0 expression counts when calculating
  correlations. This can reduce false positive correlation calculations
  when receptors have high dropout rates.

- tf_selection_method:

  Selection of which method to target transcription factors. Options:

  - "clusters": TFs are selected based on differential activity by
    cluster using Wilcoxon rank-sum test.

  - "variable": TFs are selected based on coefficient of variation (also
    known as relative standard deviation or normalized root-mean-square)
    across all cells in dataset.

  - "all": All TFs provided in the `features` matrix are included in
    downstream analysis.

- tf_variance_quantile:

  Quantile of coefficient of variation used to threshold features to
  take if `tf_selection_method` is `variable`. Features with coefficient
  of variation ranked above quantile are kept (default of 0.5 keeps the
  most variable half of features and higher numbers keep fewer
  features). Ignored if tf_selection_method is not 'variable'.

## Value

A
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object

## See also

[`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_rl_map_cellphonedb.md)
for creating receptor-ligand maps,
[`create_regulon_list_scenic()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_regulon_list_scenic.md)
for creating regulon lists from pySCENIC output, and
[`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
for building signaling networks from
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects created by this function.

## Examples

``` r

data(CellPhoneDB)
data(SCENIC)
data(PBMC)

pbmc_dom_tiny <- create_domino(
 rl_map = CellPhoneDB$rl_map_tiny, features = SCENIC$auc_tiny,
 counts = PBMC$count_tiny, z_scores = PBMC$zscore_tiny,
 clusters = PBMC$clusters_tiny, tf_targets = SCENIC$regulon_list_tiny,
 use_complexes = TRUE, remove_rec_dropout = FALSE,
 verbose = FALSE
 )
```
