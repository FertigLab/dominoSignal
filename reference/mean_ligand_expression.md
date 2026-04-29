# Calculate mean ligand expression as a data frame for plotting in circos plot

Creates a data frame of mean ligand expression for use in plotting a
circos plot of ligand expression and saving tables of mean expression.

## Usage

``` r
mean_ligand_expression(x, ligands, cell_ident, cell_barcodes, destination)
```

## Arguments

- x:

  Gene by cell expression matrix

- ligands:

  Character vector of ligand genes to be quantified

- cell_ident:

  Vector of cell type (identity) names for which to calculate mean
  ligand gene expression

- cell_barcodes:

  Vector of cell barcodes (colnames of x) belonging to cell_ident to
  calculate mean expression across

- destination:

  Name of the receptor with which each ligand interacts

## Value

A data frame of ligand expression targeting the specified receptor

## Examples

``` r
example(build_domino, echo = FALSE)
counts <- dom_counts(pbmc_dom_built_tiny)
mean_exp <- mean_ligand_expression(counts,
 ligands = c("PTPRC", "FASLG"), cell_ident = "CD14_monocyte",
 cell_barcodes = colnames(counts), destination = "FAS")
```
