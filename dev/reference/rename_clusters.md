# Renames clusters in a domino object

This function renames the clusters used to build a domino object

## Usage

``` r
rename_clusters(dom, clust_conv, warning = FALSE)
```

## Arguments

- dom:

  a domino object to rename clusters in

- clust_conv:

  named vector of conversions from old to new clusters. Values are taken
  as new clusters IDs and names as old cluster IDs.

- warning:

  logical. If TRUE, will warn if a cluster is not found in the
  conversion table. Default is FALSE.

## Value

A domino object with clusters renamed in all applicable slots.

## Examples

``` r
data(DominoObjects)
dom <- DominoObjects$built_dom_tiny
new_clust <- c(
    "CD8_T_cell" = "CD8+ T Cells",
    "CD14_monocyte" = "CD14+ Monocytes", "B_cell" = "B Cells"
)
pbmc_dom_built_tiny <- rename_clusters(dom, new_clust)
```
