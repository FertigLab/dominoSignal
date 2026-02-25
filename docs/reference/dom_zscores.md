# Access z-scores

A function to pull z-scored expression from a domino object

## Usage

``` r
dom_zscores(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

## Value

A matrix containing the z-scored gene expression values for each gene
(row) by cell (column)

## Examples

``` r
example(build_domino, echo = FALSE)
zscores <- dom_zscores(pbmc_dom_built_tiny)
```
