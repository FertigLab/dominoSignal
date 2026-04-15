# Calculate a signaling network for a domino object

This function calculates a signaling network. It requires a domino
object preprocessed from create_domino and returns a domino object
prepared for plotting with the various plotting functions in this
package.

## Usage

``` r
build_domino(
  dom,
  max_tf_per_clust = 5,
  min_tf_pval = 0.01,
  max_rec_per_tf = 5,
  rec_tf_cor_threshold = 0.15,
  min_rec_percentage = 0.1
)
```

## Arguments

- dom:

  Domino object from
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md).

- max_tf_per_clust:

  Maximum number of transcription factors called active in a cluster.

- min_tf_pval:

  Maximum p-value from differential feature score test to call a
  transcription factor active in a cluster, serving as a significance
  threshold.

- max_rec_per_tf:

  Maximum number of receptors to link to each transcription factor.

- rec_tf_cor_threshold:

  Minimum Spearman correlation used to consider a receptor linked with a
  transcription factor. Increasing this will decrease the number of
  receptors linked to each transcription factor.

- min_rec_percentage:

  Minimum percentage of cells in cluster expressing a receptor for the
  receptor to be linked to transcription factors in that cluster.

## Value

A domino object with a signaling network built

## See also

[`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
to create a domino object

## Examples

``` r
example(create_domino, echo = FALSE)

#a relaxed example
pbmc_dom_built_tiny <- build_domino(
 dom = pbmc_dom_tiny, min_tf_pval = .05, max_tf_per_clust = Inf,
 max_rec_per_tf = Inf, rec_tf_cor_threshold = .1, min_rec_percentage = 0.01
)
```
