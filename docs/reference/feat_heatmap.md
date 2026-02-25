# Create a heatmap of features organized by cluster

Creates a heatmap of transcription factor activation scores by cells
grouped by cluster.

## Usage

``` r
feat_heatmap(
  dom,
  feats = NULL,
  bool = FALSE,
  bool_thresh = 0.2,
  title = TRUE,
  norm = FALSE,
  cols = NULL,
  ann_cols = TRUE,
  min_thresh = NULL,
  max_thresh = NULL,
  ...
)
```

## Arguments

- dom:

  Domino object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md))

- feats:

  Either a vector of features to include in the heatmap or 'all' for all
  features. If left NULL then the features selected for the signaling
  network will be shown.

- bool:

  Boolean indicating whether the heatmap should be continuous or
  boolean. If boolean then bool_thresh will be used to determine how to
  define activity as positive or negative.

- bool_thresh:

  Numeric indicating the threshold separating 'on' or 'off' for feature
  activity if making a boolean heatmap.

- title:

  Either a string to use as the title or a boolean describing whether to
  include a title. In order to pass the 'main' parameter to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  you must set title to FALSE.

- norm:

  Boolean indicating whether or not to normalize the transcrption
  factors to their max value.

- cols:

  Named vector of colors to annotate cells by cluster color. Values are
  taken as colors and names as cluster. If left as NULL then default
  ggplot colors will be generated.

- ann_cols:

  Boolean indicating whether to include cell cluster as a column
  annotation. Colors can be defined with cols. If FALSE then custom
  annotations can be passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

- min_thresh:

  Minimum threshold for color scaling if not a boolean heatmap

- max_thresh:

  Maximum threshold for color scaling if not a boolean heatmap

- ...:

  Other parameters to pass to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  . Note that to use the 'main' parameter of
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  you must set title = FALSE and to use 'annCol' or 'annColors' ann_cols
  must be FALSE.

## Value

A heatmap rendered to the active graphics device

## Examples

``` r
#basic usage
example(build_domino, echo = FALSE)
feat_heatmap(pbmc_dom_built_tiny)

#using thresholds
feat_heatmap(
 pbmc_dom_built_tiny, min_thresh = 0.1, 
 max_thresh = 0.6, norm = TRUE, bool = FALSE)
#> Warning: You are using norm with min_thresh and max_thresh. Note that values will be thresholded AFTER normalization.

```
