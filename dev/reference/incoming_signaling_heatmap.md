# Create a cluster incoming signaling heatmap

Creates a heatmap of a cluster incoming signaling matrix. Each cluster
has a list of ligands capable of activating its enriched transcription
factors. The function creates a heatmap of cluster average expression
for all of those ligands. A list of all cluster incoming signaling
matrices can be found in the cl_signaling_matrices slot of a domino
option as an alternative to this plotting function.

## Usage

``` r
incoming_signaling_heatmap(
  dom,
  rec_clust,
  clusts = NULL,
  min_thresh = -Inf,
  max_thresh = Inf,
  scale = "none",
  normalize = "none",
  title = TRUE,
  ...
)
```

## Arguments

- dom:

  Domino object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

- rec_clust:

  Which cluster to select as the receptor. Must match naming of clusters
  in the domino object.

- clusts:

  Vector of clusters to be included. If NULL then all clusters are used.

- min_thresh:

  Minimum signaling threshold for plotting. Defaults to -Inf for no
  threshold.

- max_thresh:

  Maximum signaling threshold for plotting. Defaults to Inf for no
  threshold.

- scale:

  How to scale the values (after thresholding). Options are 'none',
  'sqrt' for square root, or 'log' for log10.

- normalize:

  Options to normalize the matrix. Accepted inputs are 'none' for no
  normalization, 'rec_norm' to normalize to the maximum value with each
  receptor cluster, or 'lig_norm' to normalize to the maximum value
  within each ligand cluster

- title:

  Either a string to use as the title or a boolean describing whether to
  include a title. In order to pass the 'main' parameter to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  you must set title to FALSE.

- ...:

  Other parameters to pass to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
  Note that to use the 'column_title' parameter of
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  you must set title = FALSE

## Value

a Heatmap rendered to the active graphics device

## See also

Other heatmaps:
[`cor_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/cor_heatmap.md),
[`feat_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/feat_heatmap.md),
[`signaling_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/signaling_heatmap.md)

## Examples

``` r
example(build_domino, echo = FALSE)
#incoming signaling of the CD8  T cells
incoming_signaling_heatmap(pbmc_dom_built_tiny, "CD8_T_cell")

```
