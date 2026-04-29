# Use biomaRt to convert genes

This function reads in a vector of genes and converts the genes to
specified symbol type

## Usage

``` r
convert_genes(
  genes,
  from = c("ENSMUSG", "ENSG", "MGI", "HGNC"),
  to = c("MGI", "HGNC"),
  host = "https://www.ensembl.org"
)
```

## Arguments

- genes:

  Vector of genes to convert.

- from:

  Format of gene input (ENSMUSG, ENSG, MGI, or HGNC)

- to:

  Format of gene output (MGI or HGNC)

- host:

  Host to connect to. Defaults to https://www.ensembl.org following the
  useMart default, but can be changed to archived hosts if useMart fails
  to connect.

## Value

A data frame with input genes as column 1 and converted genes as column
2
