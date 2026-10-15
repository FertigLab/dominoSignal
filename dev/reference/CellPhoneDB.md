# CellPhoneDB subset

A list of four subsets of CellPhoneDB v4.0.0 data and resulting domino
rl_map.

## Usage

``` r
data("CellPhoneDB")
```

## Format

A list of:

- genes_tiny:

  A subset of the CellPhoneDB `gene_input.csv` file.

- proteins_tiny:

  A subset of the CellPhoneDB `protein_input.csv` file.

- complexes_tiny:

  A subset of the CellPhoneDB `complex_input.csv` file.

- interactions_tiny:

  A subset of the CellPhoneDB `interaction_input.csv` file.

- rl_map_tiny:

  The subset formatted as a receptor-ligand map by
  [`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_rl_map_cellphonedb.md).

## Source

<https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz>;
subset with `data-raw/tiny_data.R`.
