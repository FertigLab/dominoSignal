# Example linkage summary

A list containing a mock linkages summary across six subjects in two
groups, and the result of running
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
on it.

## Usage

``` r
data("LinkageSummary")
```

## Format

A list of:

- linkage_sum_tiny:

  A
  [`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  object created with mock data.

- linkage_diff_tiny:

  Results from running
  [`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
  for receptors (`linkage = "rec"`) in cluster `"C1"`, grouped by
  `"group"`.

## Source

Generated with `data-raw/tiny_data.R`.
