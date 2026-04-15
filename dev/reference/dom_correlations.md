# Access correlations

A function to pull receptor-transcription factor correlations from a
domino object

## Usage

``` r
dom_correlations(dom, type = "rl")
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)

- type:

  either "rl" or "complex", to select between the receptor-ligand or
  complex correlation matrix

## Value

A matrix containing the correlation values for each receptor (row) by
transcription factor (column)

## See also

Other access:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_info.md),
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_linkages.md),
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_network_items.md),
[`dom_signaling()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_signaling.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_zscores.md)

## Examples

``` r
example(build_domino, echo = FALSE)
cor_matrix <- dom_correlations(pbmc_dom_built_tiny, "rl")
```
