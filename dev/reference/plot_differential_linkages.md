# Plot differential linkages among domino results ranked by a comparative statistic

Plot differential linkages among domino results ranked by a comparative
statistic

## Usage

``` r
plot_differential_linkages(
  differential_linkages,
  test_statistic,
  stat_range = c(0, 1),
  stat_ranking = c("ascending", "descending"),
  group_palette = NULL,
  gradient = c(minimum = "red", maximum = "gray90")
)
```

## Arguments

- differential_linkages:

  a data frame output from the
  [`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
  function

- test_statistic:

  column name of differential_linkages where the test statistic used for
  ranking linkages is stored (p.value, p.adj, odds.ratio)

- stat_range:

  a two value vector of the minimum and maximum values of test_statistic
  for plotting linkage features

- stat_ranking:

  whether to rank the features by ascending or descending values of
  test_statistic

- group_palette:

  a named vector of colors to use for each group being compared

- gradient:

  a named two value vector of colors to use for the minimum and maximum
  values of test_statistic

## Value

A heatmap-class object of features ranked by test_statistic annotated
with the proportion of subjects that showed active linkage of the
features.

## See also

Other plotting:
[`circos_ligand_receptor()`](https://FertigLab.github.io/dominoSignal/dev/reference/circos_ligand_receptor.md),
[`cor_scatter()`](https://FertigLab.github.io/dominoSignal/dev/reference/cor_scatter.md)

Differential Signaling Workflow:
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md),
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md),
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)

## Examples

``` r
data(LinkageSummary)
plot_differential_linkages(
    differential_linkages = LinkageSummary$linkage_diff_tiny,
    test_statistic = "p.value",
    stat_range = c(0, 1),
    stat_ranking = "ascending"
)

```
