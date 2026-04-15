# Access counts

A function to pull gene expression from a domino object

## Usage

``` r
dom_counts(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)

## Value

A matrix containing the gene expression values for each gene (row) by
cell (column)

## See also

Other access:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md),
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
counts <- dom_counts(pbmc_dom_built_tiny)
```
