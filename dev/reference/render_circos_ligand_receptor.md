# Render Circos Ligand Receptor Plot

Renders a circos plot from the output of
[`obtain_circos_expression()`](https://FertigLab.github.io/dominoSignal/dev/reference/obtain_circos_expression.md)
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
  [`obtain_circos_expression()`](https://FertigLab.github.io/dominoSignal/dev/reference/obtain_circos_expression.md)

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

## Examples

``` r
example(build_domino, echo = FALSE)
#basic usage
circos_df <- obtain_circos_expression(pbmc_dom_built_tiny, receptor = "CXCR3", ligands = "CCL20")
render_circos_ligand_receptor(signaling_df = circos_df, receptor = "CXCR3")
#> There are more than one numeric columns in the data frame. Take the
#> first two numeric columns and draw the link ends with unequal width.
#> 
#> Type `circos.par$message = FALSE` to suppress the message.
```
