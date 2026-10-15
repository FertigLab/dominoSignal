# Count occurrences of linkages across multiple domino results from a linkage summary

Count occurrences of linkages across multiple domino results from a
linkage summary

## Usage

``` r
count_linkage(
  link_summary,
  cluster,
  group.by = NULL,
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

  the name of the cell cluster being compared across multiple domino
  results

- group.by:

  the name of the column in `linkage_summary@subject_meta` by which to
  group subjects for counting. If NULL, only total counts of linkages
  for linkages in the cluster across all subjects is given.

- linkage:

  a stored linkage from the
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object. Can compare any of 'tfs', 'rec', 'incoming_lig', 'tfs_rec', or
  'rec_lig'

- subject_names:

  a vector of subject_names from the linkage_summary to restrict
  counting to. If NULL, all subjects in `link_summary` are included.

## Value

A data frame with columns for the unique linkage features and the counts
of how many times the linkage occurred across the compared domino
results. If group.by is used, counts of the linkages are also provided
as columns named by the unique values of the group.by variable.

## See also

[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)
to create
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
objects,
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
to compare links across groups

Differential Signaling Workflow:
[`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md),
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md),
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)

## Examples

``` r
data(LinkageSummary)
count_linkage(
  link_summary = LinkageSummary$linkage_sum_tiny, cluster = "C1", 
  group.by = "group", linkage = "rec")
#>   feature total_count G1 G2
#> 1      R1           3  3  0
#> 2      R2           4  3  1
#> 3      R3           3  1  2
#> 4      R4           2  1  1
count_linkage(
  link_summary = LinkageSummary$linkage_sum_tiny, cluster = "C1", 
  group.by = "group", linkage = "tfs_rec", subject_names = c("P1", "P2", "P4", "P5"))
#>     feature total_count G1 G2
#> 1 TF1 <- R1           2  2  0
#> 2 TF2 <- R2           3  2  1
#> 3 TF3 <- R3           3  1  2
#> 4 TF4 <- R4           2  1  1
```
