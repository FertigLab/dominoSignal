# Access clusters

A function to pull cluster information from a domino object

## Usage

``` r
dom_clusters(dom, labels = FALSE)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

- labels:

  a boolean for whether to return the cluster labels for each cell or
  the clusters used for inferring communication

## Value

A vector containing either the names of the clusters used OR factors of
the cluster label for each individual cell

## Examples

``` r
example(build_domino, echo = FALSE)
cluster_names <- dom_clusters(pbmc_dom_built_tiny)
cell_cluster_label <- dom_clusters(pbmc_dom_built_tiny, labels = TRUE)
```
