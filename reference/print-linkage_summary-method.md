# Print linkage summary object

Prints a description of a linkage summary object

## Usage

``` r
# S4 method for class 'linkage_summary'
print(x, ...)
```

## Arguments

- x:

  A linkage summary object

- ...:

  Additional arguments to be passed to other methods

## Value

A printed description of the number of subjects, groups, and clusters in
the linkage summary object

## Examples

``` r
data(LinkageSummary)
print(LinkageSummary$linkage_sum_tiny)
#> A linkage summary object of 6 subjects with 2 metadata annotations and linkages between 2 clusters.
```
