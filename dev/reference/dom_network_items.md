# Access all features, receptors, or ligands present in a signaling network.

This function collates all of the features, receptors, or ligands found
in a signaling network anywhere in a list of clusters. This can be
useful for comparing signaling networks across two separate conditions.
In order to run this
[`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
must be run on the object previously.

## Usage

``` r
dom_network_items(dom, clusters = NULL, return = NULL)
```

## Arguments

- dom:

  a domino object containing a signaling network (i.e.
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
  was run)

- clusters:

  vector indicating clusters to collate network items from. If left as
  NULL then all clusters will be included.

- return:

  string indicating whether to collate "features", "receptors", or
  "ligands". If "all" then a list of all three will be returned.

## Value

A vector containing all features, receptors, or ligands in the data set
or a list containing all three.

## See also

Other access:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_info.md),
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_linkages.md),
[`dom_signaling()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_signaling.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_zscores.md)

## Examples

``` r
example(build_domino, echo = FALSE)
monocyte_receptors <- dom_network_items(pbmc_dom_built_tiny, "CD14_monocyte", "receptors")
all_tfs <- dom_network_items(pbmc_dom_built_tiny, return = "features")
```
