# SCENIC example subset

A subset of pySCENIC output for the PBMC example data.

## Usage

``` r
data("SCENIC")
```

## Format

A list of:

- auc_tiny:

  Matrix of AUC scores from pySCENIC for a subset of TFs and cells.

- regulons_tiny:

  Data frame subset of the pySCENIC regulon table.

- regulon_list_tiny:

  The pySCENIC regulons formatted as a named list by
  [`create_regulon_list_scenic()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_regulon_list_scenic.md).

## Source

<https://doi.org/10.5281/zenodo.10124865> (files `auc_pbmc_3k.csv`,
`regulons_pbmc_3k.csv`); formatted with `data-raw/tiny_data.R`
