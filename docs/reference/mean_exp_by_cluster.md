# Get average expression for a set of genes over cluster(s)

Get average expression for a set of genes over cluster(s)

## Usage

``` r
mean_exp_by_cluster(dom, clusts, genes)
```

## Arguments

- dom:

  A domino object

- clusts:

  Cluster(s) for which we want to get average expression

- genes:

  The genes for which we want to get average expression

## Value

A dataframe of genes x clusters, values are z-scores averaged over the
clusters
