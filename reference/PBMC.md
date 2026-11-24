# PBMC scRNAseq data subset

A subset of the preprocessed 10x Genomics PBMC 3k scRNAseq data.

## Usage

``` r
data("PBMC")
```

## Format

A list of:

- count_tiny:

  Subset of the sparse counts matrix (`dgCMatrix`), genes x cells.

- zscore_tiny:

  Subset of z-scored expression matrix, genes x cells.

- clusters_tiny:

  Factor of cell-type labels named by cell barcode.

## Source

<https://doi.org/10.5281/zenodo.10124865> (file `pbmc3k_sce.rds`);
formatted with `data-raw/tiny_data.R`.
