# Obtain Circos Expression

Pull expression data from a domino object and format for plotting as a
receptor-oriented circos plot.

## Usage

``` r
obtain_circos_expression(
  dom,
  receptor,
  ligands,
  ligand_expression_threshold = 0.01,
  cell_idents = NULL
)
```

## Arguments

- dom:

  Domino object that has undergone network building with build_domino()

- receptor:

  Name of a receptor active in at least one cell type in the domino
  object

- ligands:

  Character vector of ligands capable of interaction with the receptor

- ligand_expression_threshold:

  Minimum mean expression value of a ligand by a cell type for a chord
  to be rendered between the cell type and the receptor

- cell_idents:

  Vector of cell types from cluster assignments in the domino object to
  be included in the plot.

## Value

a data frame where each row describes plotting parameters of
ligand-receptor interactions to pass to render_circos_ligand_receptor()
