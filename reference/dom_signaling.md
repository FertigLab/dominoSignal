# Access signaling

A function to pull signaling matrices from a domino object

## Usage

``` r
dom_signaling(dom, cluster = NULL)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

- cluster:

  either NULL to indicate global signaling or a specific cluster for
  which a signaling matrix is desired

## Value

A data frame containing the signaling score through each ligand (row) by
each cluster (column) OR a data frame containing the global summed
signaling scores between receptors (rows) and ligands (columns) of each
cluster

## See also

Accessor Functions:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/reference/dom_info.md),
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/reference/dom_linkages.md),
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/reference/dom_network_items.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/reference/dom_zscores.md)

## Examples

``` r
example(build_domino, echo = FALSE)
monocyte_signaling <- dom_signaling(pbmc_dom_built_tiny, cluster = "CD14_monocyte")
```
