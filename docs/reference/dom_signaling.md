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

## Examples

``` r
example(build_domino, echo = FALSE)
monocyte_signaling <- dom_signaling(pbmc_dom_built_tiny, cluster = "CD14_monocyte")
```
