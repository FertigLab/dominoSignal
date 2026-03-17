# Access correlations

A function to pull receptor-transcription factor correlations from a
domino object

## Usage

``` r
dom_correlations(dom, type = "rl")
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

- type:

  either "rl" or "complex", to select between the receptor-ligand or
  complex correlation matrix

## Value

A matrix containing the correlation values for each receptor (row) by
transcription factor (column)

## Examples

``` r
example(build_domino, echo = FALSE)
cor_matrix <- dom_correlations(pbmc_dom_built_tiny, "rl")
```
