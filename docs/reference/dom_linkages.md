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
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)

- link_type:

  one value (out of "complexes", "receptor-ligand", "tf-target",
  "tf-receptor", "receptor", "incoming-ligand") used to select the
  desired type of linkage

- by_cluster:

  a boolean to indicate whether the linkages should be returned overall
  or by cluster

## Value

A list containing linkages between some combination of receptors,
ligands, transcription factors, and clusters

## Examples

``` r
example(build_domino, echo = FALSE)
complexes <- dom_linkages(pbmc_dom_built_tiny, "complexes")
tf_rec_by_cluster <- dom_linkages(pbmc_dom_built_tiny, "tf-receptor", TRUE)
```
