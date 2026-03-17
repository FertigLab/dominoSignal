# Count occurrences of linkages across multiple domino results from a linkage summary

Count occurrences of linkages across multiple domino results from a
linkage summary

## Usage

``` r
count_linkage(
  linkage_summary,
  cluster,
  group.by = NULL,
  linkage = "rec_lig",
  subject_names = NULL
)
```

## Arguments

- linkage_summary:

  a
  [`linkage_summary()`](https://FertigLab.github.io/dominoSignal/reference/linkage_summary-class.md)
  object

- cluster:

  the name of the cell cluster being compared across multiple domino
  results

- group.by:

  the name of the column in `linkage_summary@subject_meta` by which to
  group subjects for counting. If NULL, only total counts of linkages
  for linkages in the cluster across all subjects is given.

- linkage:

  a stored linkage from the domino object. Can compare any of 'tfs',
  'rec', 'incoming_lig', 'tfs_rec', or 'rec_lig'

- subject_names:

  a vector of subject_names from the linkage_summary to be compared. If
  NULL, all subject_names in the linkage summary are included in
  counting.

## Value

A data frame with columns for the unique linkage features and the counts
of how many times the linkage occured across the compared domino
results. If group.by is used, counts of the linkages are also provided
as columns named by the unique values of the group.by variable.

## Examples

``` r
count_linkage(
  linkage_summary = mock_linkage_summary(), cluster = "C1", 
  group.by = "group", linkage = "rec")
#>   feature total_count G1 G2
#> 1      R1           3  3  0
#> 2      R2           4  3  1
#> 3      R3           3  1  2
#> 4      R4           2  1  1
```
