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

## Examples

``` r
example(build_domino, echo = FALSE)
#basic usage
obtain_circos_expression(pbmc_dom_built_tiny, receptor = "CXCR3", ligands = "CCL20")
#>                origin destination mean.expression        sender ligand receptor
#> 1        B_cell-CCL20       CXCR3       0.0000000        B_cell  CCL20    CXCR3
#> 2 CD14_monocyte-CCL20       CXCR3       0.2095956 CD14_monocyte  CCL20    CXCR3
#> 3    CD8_T_cell-CCL20       CXCR3       0.0000000    CD8_T_cell  CCL20    CXCR3
#>   scaled.mean.expression ligand.arc receptor.arc
#> 1                      0          1     1.333333
#> 2                      1          1     1.333333
#> 3                      0          1     1.333333
```
