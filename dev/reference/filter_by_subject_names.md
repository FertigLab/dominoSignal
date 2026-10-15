# Restrict a linkage_summary object to a set of subject names

Restrict a linkage_summary object to a set of subject names

## Usage

``` r
filter_by_subject_names(link_summary, subject_names)
```

## Arguments

- link_summary:

  a
  [`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  object

- subject_names:

  a vector of subject_names to restrict `link_summary` to, or NULL to
  leave `link_summary` unchanged

## Value

`link_summary`, subset to the given `subject_names`, or unchanged if
`subject_names` is NULL
