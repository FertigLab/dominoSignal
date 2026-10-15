# Create a correlation plot between TF and receptor

Create a correlation plot between transcription factor activity scores
and receptor expression

## Usage

``` r
cor_scatter(dom, tf, rec, remove_rec_dropout = NULL, ...)
```

## Arguments

- dom:

  A
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
  has been run)

- tf:

  Target TF for plotting AUC score

- rec:

  Target receptor for plotting expression

- remove_rec_dropout:

  Whether to remove cells with zero expression for plot. Default (NULL)
  will use the value of `remove_rec_dropout` used with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
  when the
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object was created. If the
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object was created before that parameter was added, it will default to
  FALSE. To override the default, set this parameter to TRUE or FALSE

- ...:

  Other parameters to pass to
  [`ggpubr::ggscatter()`](https://rpkgs.datanovia.com/ggpubr/reference/ggscatter.html).

## Value

A ggplot scatter plot rendered in the active graphics device

## See also

Other plotting:
[`circos_ligand_receptor()`](https://FertigLab.github.io/dominoSignal/dev/reference/circos_ligand_receptor.md),
[`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md)

## Examples

``` r
data(DominoObjects)
dom <- DominoObjects$built_dom_tiny
cor_scatter(dom, "FLI1","CXCR3")

```
