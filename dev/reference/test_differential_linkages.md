# Statistical test for differential linkages across multiple domino results

Statistical test for differential linkages across multiple domino
results

## Usage

``` r
test_differential_linkages(
  link_summary,
  cluster,
  group.by,
  linkage = "rec_lig",
  subject_names = NULL
)
```

## Arguments

- link_summary:

  a
  [`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  object

- cluster:

  the name of the cell cluster being compared across multiple
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  results

- group.by:

  the name of the column in `linkage_summary@subject_meta` by which to
  group subjects for counting.

- linkage:

  a stored linkage from the
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object. Can compare any of 'tfs', 'rec', 'incoming_lig', 'tfs_rec', or
  'rec_lig'

- subject_names:

  a vector of subject_names from the linkage_summary to restrict testing
  to. If NULL, all subjects in `link_summary` are included.

## Value

A data frame of results from the test of the differential linkages. Rows
correspond to each linkage tested. Columns correspond to:

- 'cluster' : the name of the cell cluster being compared

- 'linkage' : the type of linkage being compared

- 'group.by' : the grouping variable

- 'test_name' : the test used for comparison (Fisher's exact test)

- 'feature' : individual linkages compared

- 'test statistics' : an odds ratio (for 2 groups), p-value, and
  fdr-adjusted p-value.

- 'total_count' : total number of subjects where the linkage is active

- 'X_count' : number of subjects in each category of group.by (X) where
  the linkage is active

- 'total_n' : number of total subjects compared

- 'X_n' : total number of subjects in each category of group.by (X)

## See also

[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)
to create
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
objects

Differential Signaling Workflow:
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md),
[`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md),
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)

## Examples

``` r
data(LinkageSummary)
tiny_differential_linkage_c1 <- test_differential_linkages(
    link_summary = LinkageSummary$linkage_sum_tiny, cluster = "C1", group.by = "group",
    linkage = "rec", subject_names = c("P1", "P2", "P4", "P5")
)
```
