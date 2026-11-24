# Access signaling

A function to pull signaling matrices from a
[`domino()`](https://FertigLab.github.io/dominoSignal/reference/domino-class.md)
object `dom_signaling()` and `dom_signalling()` are synonyms.

## Usage

``` r
dom_signaling(dom, cluster = NULL)

dom_signalling(dom, cluster = NULL)
```

## Arguments

- dom:

  a
  [`domino()`](https://FertigLab.github.io/dominoSignal/reference/domino-class.md)
  object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

- cluster:

  either NULL to indicate global signaling or a specific cluster for
  which a signaling matrix is desired

## Value

A data frame containing the signaling score through each ligand (row) by
each cluster (column) OR . a data frame containing the global summed
signaling scores between receptors (rows) and ligands (columns) of each
cluster

## See also

Other accessors:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/reference/dom_info.md),
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/reference/dom_linkages.md),
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/reference/dom_network_items.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/reference/dom_zscores.md)

## Examples

``` r
data(DominoObjects)
dom <- DominoObjects$built_dom_tiny
monocyte_signaling <- dom_signaling(dom, cluster = "CD14_monocyte")
```
