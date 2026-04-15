# Access linkages

A function to pull linkages from a domino object

## Usage

``` r
dom_linkages(
  dom,
  link_type = c("complexes", "receptor-ligand", "tf-target", "tf-receptor", "receptor",
    "incoming-ligand"),
  by_cluster = FALSE
)
```

## Arguments

- dom:

  a domino object that has been created with
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)

- link_type:

  one value (out of "complexes", "receptor-ligand", "tf-target",
  "tf-receptor", "receptor", "incoming-ligand") used to select the
  desired type of linkage. Note that "receptor" and "incoming-ligand"
  are only available when by_cluster is set to TRUE.

- by_cluster:

  a boolean to indicate whether the linkages should be returned overall
  or by cluster

## Value

A list containing linkages between some combination of receptors,
ligands, transcription factors, and clusters

## See also

Other access:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md),
[`dom_database()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_database.md),
[`dom_de()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_de.md),
[`dom_info()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_info.md),
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_network_items.md),
[`dom_signaling()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_signaling.md),
[`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_tf_activation.md),
[`dom_zscores()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_zscores.md)

## Examples

``` r
example(build_domino, echo = FALSE)
complexes <- dom_linkages(pbmc_dom_built_tiny, "complexes")
tf_rec_by_cluster <- dom_linkages(pbmc_dom_built_tiny, "tf-receptor", TRUE)
```
