# Print domino object

Prints a summary of a domino object

## Usage

``` r
# S4 method for class 'domino'
print(x, ...)
```

## Arguments

- x:

  A domino object

- ...:

  Additional arguments to be passed to other methods

## Value

A printed description of the number of cells and clusters in the domino
object

## Examples

``` r
example(build_domino, echo = FALSE)
print(pbmc_dom_built_tiny)
#> A domino object of 360 cells
#>                 Contains signaling between 3 clusters
#>                 Built with a maximum of Inf TFs per cluster
#>                 and a maximum of Inf receptors per TF
```
