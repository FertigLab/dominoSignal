# Get ligand-receptor signaling information

Get ligand-receptor signaling information

## Usage

``` r
get_signaling_info(dom, rec_clusters, cl_ligands_sub, exp_type)
```

## Arguments

- dom:

  A built domino object (as output by
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

- rec_clusters:

  Character or factor vector of cluster names for which to compute
  incoming receptor/TF signals

- cl_ligands_sub:

  Data frame with columns 'ligand', 'cluster', 'mean_counts'; using
  [`reshape2::melt()`](https://rdrr.io/pkg/reshape2/man/melt.html) on
  the output of
  [`get_ligand_expression()`](https://FertigLab.github.io/dominoSignal/dev/reference/get_ligand_expression.md)
  is a convenient way to get this

- exp_type:

  Character of length 1: either "counts" or "z_scores" to specify
  expression type

## Value

Data frame with columns: ligand, receptor, transcription_factor,
ligand_exp, rec_exp, tf_auc, sending_cl, receiving_cl. Each row
represents a ligand-receptor-TF triplet with receiver cluster, and
corresponding expression values. Returns empty data frame if no
interactions found.
