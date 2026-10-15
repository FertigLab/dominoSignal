# Subset a linkage_summary object

Subsets a linkage summary object by subject names or metadata

## Usage

``` r
# S4 method for class 'linkage_summary'
subset(x, subset)
```

## Arguments

- x:

  A linkage_summary object

- subset:

  A logical expression referencing `subject_names` or columns of
  `subject_meta`

## Value

A linkage_summary object containing only the subjects that match the
specified criteria.

## Examples

``` r
data(LinkageSummary)
links <- LinkageSummary$linkage_sum_tiny
subset(links, subset = !subject_names %in% c("P1", "P2"))
#> A linkage summary object of 4 subjects with 2 metadata annotations and linkages between 2 clusters.
subset(links, subset = group == "G1")
#> A linkage summary object of 3 subjects with 2 metadata annotations and linkages between 2 clusters.
subset(links, subset = subject_names %in% c("P1", "P2", "P4") | group == "G2")
#> A linkage summary object of 5 subjects with 2 metadata annotations and linkages between 2 clusters.
```
