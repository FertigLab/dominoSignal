# Show linkage_summary object information

Shows content overview of linkage_summary object

## Usage

``` r
# S4 method for class 'linkage_summary'
show(object)
```

## Arguments

- object:

  A linkage_summary object

## Value

A printed description of the number of subjects, groups, and clusters in
the linkage summary object

## Examples

``` r
link_sum <- mock_linkage_summary()
show(link_sum)
#> A linkage summary object of 6 subjects with 2 metadata annotations and linkages between 2 clusters.
```
