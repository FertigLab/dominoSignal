# Access database

A function to pull database information from a domino object

## Usage

``` r
dom_database(dom, name_only = TRUE)
```

## Arguments

- dom:

  a domino object that has been created

- name_only:

  a boolean for whether to return only the name of the database used or
  the entire database that is stored. Default TRUE.

## Value

A vector of unique databases used in building the domino object OR a
data frame that includes the database information used in the domino
object creation

## See also

Other access:
[`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md),
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md),
[`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md),
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
database_name <- dom_database(pbmc_dom_built_tiny)
full_database <- dom_database(pbmc_dom_built_tiny, name_only = FALSE)
```
