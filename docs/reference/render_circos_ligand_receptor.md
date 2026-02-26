# Render Circos Ligand Receptor Plot

Renders a circos plot from the output of
[`obtain_circos_expression()`](https://FertigLab.github.io/dominoSignal/reference/obtain_circos_expression.md)
to the active graphics device

## Usage

``` r
render_circos_ligand_receptor(
  signaling_df,
  receptor,
  cell_colors = NULL,
  ligand_expression_threshold = 0.01
)
```

## Arguments

- signaling_df:

  Data frame output from
  [`obtain_circos_expression()`](https://FertigLab.github.io/dominoSignal/reference/obtain_circos_expression.md)

- receptor:

  Name of a receptor active in at least one cell type in the domino
  object

- cell_colors:

  Named vector of color names or hex codes where names correspond to the
  plotted cell types and the color values

- ligand_expression_threshold:

  Minimum mean expression value of a ligand by a cell type for a chord
  to be rendered between the cell type and the receptor

## Value

a circlize plot is rendered to the active graphics device
