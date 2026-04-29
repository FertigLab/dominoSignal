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
link_sum <- mock_linkage_summary()
print(link_sum)
#> A linkage summary object of 6 subjects with 2 metadata annotations and linkages between 2 clusters.
```
