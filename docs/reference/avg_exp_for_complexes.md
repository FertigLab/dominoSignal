# Get average expression for complexes

Get average expression for complexes

## Usage

``` r
avg_exp_for_complexes(exp_mat, complexes_list)
```

## Arguments

- exp_mat:

  A matrix(or dataframe) of genes x clusters, values are z-scores
  averaged over the clusters

- complexes_list:

  A list similar to dom@linkages\$complexes

## Value

A list containing average expression for any complexes
