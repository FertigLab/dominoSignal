# Get ligand expression matrix for outgoing clusters

Get ligand expression matrix for outgoing clusters

## Usage

``` r
get_ligand_expression(dom, send_clusters, lig_genes, complexes, exp_type)
```

## Arguments

- dom:

  A built domino object (as output by
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

- send_clusters:

  Character or factor vector of cluster names for which to compute
  outgoing ligand signals

- lig_genes:

  Character vector of ligand gene names; should intersect with
  expression matrix rownames

- complexes:

  Named list where names are complex identifiers and values are
  character vectors of component genes. Used to aggregate expression for
  multi-gene complexes. Empty list is acceptable. The output of
  [`get_resolved_ligands()`](https://FertigLab.github.io/dominoSignal/dev/reference/get_resolved_ligands.md)
  can be used here.

- exp_type:

  Character of length 1: either "counts" or "z_scores" to specify
  expression type

## Value

Matrix with ligands/complexes as rows and send_clusters as columns;
entries are mean expression. Missing ligands and empty clusters yield
NA; row and column names are preserved.
