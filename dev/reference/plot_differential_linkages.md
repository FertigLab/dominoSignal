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
  group_palette = NULL
)
```

## Arguments

- differential_linkages:

  a data frame output from the
  [`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
  function

- test_statistic:

  column name of differential_linkages where the test statistic used for
  ranking linkages is stored (ex. 'p.value')

- stat_range:

  a two value vector of the minimum and maximum values of test_statistic
  for plotting linkage features

- stat_ranking:

  'ascending' (lowest value of test statisic is colored red and plotted
  at the top) or 'descending' (highest value of test statistic is
  colored red and plotted at the top).

- group_palette:

  a named vector of colors to use for each group being compared

## Value

A heatmap-class object of features ranked by test_statistic annotated
with the proportion of subjects that showed active linkage of the
features.

## See also

Other misc_plotting:
[`circos_ligand_receptor()`](https://FertigLab.github.io/dominoSignal/dev/reference/circos_ligand_receptor.md),
[`cor_scatter()`](https://FertigLab.github.io/dominoSignal/dev/reference/cor_scatter.md)

Other differentials:
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md),
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md),
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)

## Examples

``` r
example(build_domino, echo = FALSE)
example(test_differential_linkages, echo = FALSE)
plot_differential_linkages(
 differential_linkages = tiny_differential_linkage_c1,
 test_statistic = "p.value",
 stat_ranking = "ascending"
)

```
