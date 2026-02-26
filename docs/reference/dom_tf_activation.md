# Access transcription factor activation

A function to pull transcription factor activation scores from a domino
object

## Usage

``` r
dom_tf_activation(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

## Value

A matrix containing the transcription factor activation scores for each
TF (row) by cell (column)

## Examples

``` r
example(build_domino, echo = FALSE)
tf_activation <- dom_tf_activation(pbmc_dom_built_tiny)
```
