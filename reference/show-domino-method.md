# Show domino object information

Shows content overview of domino object

## Usage

``` r
# S4 method for class 'domino'
show(object)
```

## Arguments

- object:

  A domino object

## Value

A printed description of cell numbers and clusters in the object

## Examples

``` r
example(build_domino, echo = FALSE)
show(pbmc_dom_built_tiny)
#> A domino object of 360 cells
#> Built with signaling between 3 clusters
```
