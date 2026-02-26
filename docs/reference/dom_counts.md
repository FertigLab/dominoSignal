# Access counts

A function to pull gene expression from a domino object

## Usage

``` r
dom_counts(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

## Value

A matrix containing the gene expression values for each gene (row) by
cell (column)

## Examples

``` r
example(build_domino, echo = FALSE)
counts <- dom_counts(pbmc_dom_built_tiny)
```
