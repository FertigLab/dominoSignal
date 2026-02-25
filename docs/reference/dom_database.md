# Access database

A function to pull database information from a domino object

## Usage

``` r
dom_database(dom, name_only = TRUE)
```

## Arguments

- dom:

  a domino object that has been created

- name_only:

  a boolean for whether to return only the name of the database used or
  the entire database that is stored. Default TRUE.

## Value

A vector of unique databases used in building the domino object OR a
data frame that includes the database information used in the domino
object creation

## Examples

``` r
example(build_domino, echo = FALSE)
database_name <- dom_database(pbmc_dom_built_tiny)
full_database <- dom_database(pbmc_dom_built_tiny, name_only = FALSE)
```
