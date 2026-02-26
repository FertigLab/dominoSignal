# Create a network heatmap

Creates a heatmap of the signaling network. Alternatively, the network
matrix can be accessed directly in the signaling slot of a domino object
using the
[`dom_signaling()`](https://FertigLab.github.io/dominoSignal/reference/dom_signaling.md)
function.

## Usage

``` r
signaling_heatmap(
  dom,
  clusts = NULL,
  min_thresh = -Inf,
  max_thresh = Inf,
  scale = "none",
  normalize = "none",
  ...
)
```

## Arguments

- dom:

  domino object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md))

- clusts:

  vector of clusters to be included. If NULL then all clusters are used.

- min_thresh:

  minimum signaling threshold for plotting. Defaults to -Inf for no
  threshold.

- max_thresh:

  maximum signaling threshold for plotting. Defaults to Inf for no
  threshold.

- scale:

  how to scale the values (after thresholding). Options are 'none',
  'sqrt' for square root, or 'log' for log10.

- normalize:

  options to normalize the matrix. Normalization is done after
  thresholding and scaling. Accepted inputs are 'none' for no
  normalization, 'rec_norm' to normalize to the maximum value with each
  receptor cluster, or 'lig_norm' to normalize to the maximum value
  within each ligand cluster

- ...:

  other parameters to pass to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Value

A heatmap rendered to the active graphics device

## Examples

``` r
example(build_domino, echo = FALSE)
#basic usage
signaling_heatmap(pbmc_dom_built_tiny)

#scale
signaling_heatmap(pbmc_dom_built_tiny, scale = "sqrt")

#normalize
signaling_heatmap(pbmc_dom_built_tiny, normalize = "rec_norm")
#> Warning: Some values are NA, replacing with 0s.

```
