# Create a correlation plot between TF and receptor

Create a correlation plot between transcription factor activation score
and receptor expression

## Usage

``` r
cor_scatter(dom, tf, rec, remove_rec_dropout = TRUE, ...)
```

## Arguments

- dom:

  Domino object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md))

- tf:

  Target TF for plottting AUC score

- rec:

  Target receptor for plotting expression

- remove_rec_dropout:

  Whether to remove cells with zero expression for plot. This should
  match the same setting as in
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md).

- ...:

  Other parameters to pass to
  [`ggpubr::ggscatter()`](https://rpkgs.datanovia.com/ggpubr/reference/ggscatter.html).

## Value

A ggplot scatter plot rendered in the active graphics device

## Examples

``` r
example(build_domino, echo = FALSE)
cor_scatter(pbmc_dom_built_tiny, "FLI1","CXCR3")

```
