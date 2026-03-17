# The domino class

The domino class contains all information necessary to calculate
receptor-ligand signaling. It contains z-scored expression, cell cluster
labels, feature values, and a referenced receptor-ligand database
formatted as a receptor-ligand map. Calculated intermediate values are
also stored.

## Value

An instance of class `domino `

## Slots

- `db_info`:

  List of data sets from ligand - receptor database

- `counts`:

  Raw count gene expression data

- `z_scores`:

  Matrix of z-scored expression data with cells as columns

- `clusters`:

  Named factor with cluster identity of each cell

- `features`:

  Matrix of features (TFs) to correlate receptor - ligand expression
  with. Cells are columns and features are rows.

- `cor`:

  Correlation matrix of receptor expression to features.

- `linkages`:

  List of lists containing info linking cluster-\>tf-\>rec-\>lig

- `clust_de`:

  Data frame containing differential expression results for features by
  cluster.

- `misc`:

  List of miscellaneous info pertaining to run parameters etc.

- `cl_signaling_matrices`:

  Incoming signaling matrix for each cluster

- `signaling`:

  Signaling matrix between all clusters.
