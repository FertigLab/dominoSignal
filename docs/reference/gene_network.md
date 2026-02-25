# Create a gene association network

Create a gene association network for genes from a given cluster. The
selected cluster acts as the receptor for the gene association network,
so only ligands, receptors, and features associated with the receptor
cluster will be included in the plot.

## Usage

``` r
gene_network(
  dom,
  clust = NULL,
  OutgoingSignalingClust = NULL,
  class_cols = c(lig = "#FF685F", rec = "#47a7ff", feat = "#39C740"),
  cols = NULL,
  lig_scale = 1,
  layout = "grid",
  ...
)
```

## Arguments

- dom:

  Domino object with network built
  ([`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md))

- clust:

  Receptor cluster to create the gene association network for. A vector
  of clusters may be provided.

- OutgoingSignalingClust:

  Vector of clusters to plot the outgoing signaling from

- class_cols:

  Named vector of colors used to color classes of vertices. Values must
  be colors and names must be classes ('rec', 'lig', and 'feat' for
  receptors, ligands, and features.).

- cols:

  Named vector of colors for individual genes. Genes not included in
  this vector will be colored according to class_cols.

- lig_scale:

  FALSE or a numeric value to scale the size of ligand vertices based on
  z-scored expression in the data set.

- layout:

  Type of layout to use. Options are 'grid', 'random', 'sphere',
  'circle', 'fr' for Fruchterman-Reingold force directed layout, and
  'kk' for Kamada Kawai for directed layout.

- ...:

  Other parameters to pass to plot() with an
  [igraph](https://r.igraph.org/) object. See
  [igraph](https://r.igraph.org/) manual for options.

## Value

An igraph plot rendered to the active graphics device

## Examples

``` r
#basic usage
example(build_domino, echo = FALSE)
gene_network(
 pbmc_dom_built_tiny, clust = "CD8_T_cell", 
 OutgoingSignalingClust = "CD14_monocyte")

```
