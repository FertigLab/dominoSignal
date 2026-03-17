# Access differential expression

A function to pull differential expression p-values from a domino object

## Usage

``` r
dom_de(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

## Value

A matrix containing the p-values for differential expression of
transcription factors (rows) in each cluster (columns)

## Examples

``` r
example(build_domino, echo = FALSE)
de_mat <- dom_de(pbmc_dom_built_tiny)
```
