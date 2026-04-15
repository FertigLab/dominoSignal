# Statistical test for differential linkages across multiple domino results

Statistical test for differential linkages across multiple domino
results

## Usage

``` r
test_differential_linkages(
  linkage_summary,
  cluster,
  group.by,
  linkage = "rec_lig",
  subject_names = NULL,
  test_name = "fishers.exact"
)
```

## Arguments

- linkage_summary:

  a
  [`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  object

- cluster:

  the name of the cell cluster being compared across multiple domino
  results

- group.by:

  the name of the column in `linkage_summary@subject_meta` by which to
  group subjects for counting.

- linkage:

  a stored linkage from the domino object. Can compare any of 'tfs',
  'rec', 'incoming_lig', 'tfs_rec', or 'rec_lig'

- subject_names:

  a vector of subject_names from the linkage_summary. NOTE: all
  subject_names in the linkage summary are included in counting.

- test_name:

  the statistical test used for comparison.

  - 'fishers.exact' : Fisher's exact test for the dependence of the
    proportion of subjects with an active linkage in the cluster on
    which group the subject belongs to in the group.by variable.
    Provides an odds ratio, p-value, and a Benjamini-Hochberg
    FDR-adjusted p-value (p.adj) for each linkage tested.

## Value

A data frame of results from the test of the differential linkages. Rows
correspond to each linkage tested. Columns correspond to:

- 'cluster' : the name of the cell cluster being compared

- 'linkage' : the type of linkage being compared

- 'group.by' : the grouping variable

- 'test_name' : the test used for comparison

- 'feature' : individual linkages compared

- 'test statistics' : test statistics provided are based on test method.
  'fishers.exact' provides a odds ratio, p-value, and fdr-adjusted
  p-value.

- 'total_count' : total number of subjects where the linkage is active

- 'X_count' : number of subjects in each category of group.by (X) where
  the linkage is active

- 'total_n' : number of total subjects compared

- 'X_n' : total number of subjects in each category of group.by (X)

## See also

Other differentials:
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md),
[`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md),
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)

## Examples

``` r
tiny_differential_linkage_c1 <- test_differential_linkages(
  linkage_summary = mock_linkage_summary(), cluster = "C1", group.by = "group",
  linkage = "rec", test_name = "fishers.exact"
)
```
