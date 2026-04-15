# Convert genes using a table

Takes a vector of gene inputs and a conversion table and returns a
converted gene table

## Usage

``` r
table_convert_genes(genes, from, to, conversion_table)
```

## Arguments

- genes:

  the genes to convert

- from:

  gene symbol type of the input (one of ENSG, ENSMUSG, HGNC, MGI)

- to:

  desired gene symbol type for the output (one of HGNC, MGI)

- conversion_table:

  a data frame with column names corresponding to gene symbol types
  (mm.ens, hs.ens, mgi, hgnc) and rows corresponding to the gene symbols
  themselves

## Value

A data frame of genes with original and corresponding converted symbols
