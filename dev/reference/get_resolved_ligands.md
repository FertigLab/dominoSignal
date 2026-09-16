# Get ligands with resolved names

Get ligands with resolved names

## Usage

``` r
get_resolved_ligands(dom)
```

## Arguments

- dom:

  A built domino object (as output by
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

## Value

A named list with two elements:

- lig_names:

  Character vector of unique ligand gene names or aliases

- complex_names:

  Named list mapping complex names to component gene vectors
