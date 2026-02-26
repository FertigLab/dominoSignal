# Check input arguments

Accepts an object and rules to check against; stops if requirements are
not met

## Usage

``` r
check_arg(
  arg,
  allow_class = NULL,
  allow_len = NULL,
  allow_range = NULL,
  allow_values = NULL,
  need_vars = c(NULL),
  need_colnames = FALSE,
  need_rownames = FALSE,
  need_names = FALSE
)
```

## Arguments

- arg:

  the argument to check

- allow_class:

  vector of allowed classes

- allow_len:

  vector of allowed lengths

- allow_range:

  range of minimum and maximum values i.e. c(1, 5)

- allow_values:

  vector of allowed values

- need_vars:

  vector of required variables

- need_colnames:

  vogical for whether colnames are required

- need_rownames:

  logical for whether rownames are required

- need_names:

  logical for whether names are required

## Value

Logical indicating whether the argument meets the requirements
