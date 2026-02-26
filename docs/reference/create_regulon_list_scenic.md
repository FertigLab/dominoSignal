# Create a list of genes in regulons inferred by SCENIC

Generates a list of transcription factors and the genes targeted by the
transcription factor as part of their regulon inferred by pySCENIC

## Usage

``` r
create_regulon_list_scenic(regulons)
```

## Arguments

- regulons:

  Data frame or file path to the table of the output of the ctx function
  from pySCENIC

## Value

A list where names are transcription factors and the stored values are
character vectors of genes in the inferred regulons

## Examples

``` r
data(SCENIC)
regulon_list_tiny <- create_regulon_list_scenic(regulons = SCENIC$regulons_tiny)
```
