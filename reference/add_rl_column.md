# Adds a column to the RL signaling data frame.

This function adds a column to the internal rl 'map' used to map all
receptor and receptor complexes to all ligand and ligand complexes.

## Usage

``` r
add_rl_column(map, map_ref, conv, new_name)
```

## Arguments

- map:

  RL signaling data frame.

- map_ref:

  Name of column to match new data to

- conv:

  Data frame matching current data in map to new data.

- new_name:

  Name of new column to be created in RL map

## Value

An updated RL signaling data frame

## Examples

``` r
example(create_rl_map_cellphonedb, echo = FALSE)
lr_name <- data.frame("abbrev" = c("L", "R"), "full" = c("Ligand", "Receptor"))
rl_map_expanded <- add_rl_column(map = rl_map_tiny, map_ref = "type_A",
conv = lr_name, new_name = "type_A_full")
```
