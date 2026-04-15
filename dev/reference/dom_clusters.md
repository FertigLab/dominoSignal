# Access clusters

A function to pull cluster information from a domino object

## Usage

``` r
dom_clusters(dom, labels = FALSE)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)

- labels:

  a boolean for whether to return the cluster labels for each cell or
  the clusters used for inferring communication

## Value

A vector containing either the names of the clusters used OR factors of
the cluster label for each individual cell

## See also

Other access:
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_info.md),
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_linkages.md),
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_network_items.md),
[`dom_signaling()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_signaling.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_zscores.md)

## Examples

``` r
example(build_domino, echo = FALSE)
cluster_names <- dom_clusters(pbmc_dom_built_tiny)
cell_cluster_label <- dom_clusters(pbmc_dom_built_tiny, labels = TRUE)
```
