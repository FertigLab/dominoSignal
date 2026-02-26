# Convert between complex names and gene names

Convert between complex names and gene names

## Usage

``` r
resolve_complexes(dom, genes)
```

## Arguments

- dom:

  A domino object

- genes:

  A vector of genes, some of which may be complexes

## Value

A list where any complexes are mapped to a vector of component genes.
The list names are set to the input gene names.
