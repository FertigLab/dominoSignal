# Access build information

A function to pull the parameters used when running
[`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md)
from a domino object

## Usage

``` r
dom_info(dom)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

## Value

A list containing booleans for whether the object has been created and
built and a list of the build parameters that were used in
[`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md)
to infer the signaling network

## Examples

``` r
example(build_domino, echo = FALSE)
build_details <- dom_info(pbmc_dom_built_tiny)
```
