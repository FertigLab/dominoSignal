# Show domino object information

Shows content overview of a
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object

## Usage

``` r
# S4 method for class 'domino'
show(object)
```

## Arguments

- object:

  A
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object

## Value

A printed description of cell numbers and clusters in the object

## Examples

``` r
data(DominoObjects)
show(DominoObjects$built_dom_tiny)
#> A domino object of 360 cells
#>                 Contains signaling between 3 clusters
#>                 Built with a maximum of Inf TFs per cluster
#>                 and a maximum of Inf receptors per TF
#> Created using 'clusters' TF selection with receptor and ligand complexes
#> Created with dominoSignal v1.7.1
#> Built with dominoSignal v1.7.1
```
