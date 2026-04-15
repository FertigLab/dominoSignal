# Package index

## Domino Classes

An object class designed to store information used and produced in
dominoSignal analysis

- [`domino-class`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  [`domino`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  : The domino class
- [`linkage_summary-class`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  [`linkage_summary`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
  : The domino linkage summary class

## Analysis Functions

Functions for performing cell-cell communication analysis with the
dominoSignal package

- [`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_rl_map_cellphonedb.md)
  : Create a receptor - ligand map from a CellPhoneDB signaling database
- [`create_regulon_list_scenic()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_regulon_list_scenic.md)
  : Create a list of genes in regulons inferred by SCENIC
- [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
  : Create a domino object and prepare it for network construction
- [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
  : Calculate a signaling network for a domino object
- [`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)
  : Summarize linkages from multiple domino objects
- [`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
  : Statistical test for differential linkages across multiple domino
  results

## Plotting Functions

Functions for visualizing dominoSignal analysis results

### Heatmaps

- [`cor_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/cor_heatmap.md)
  : Create a heatmap of correlation between receptors and transcription
  factors
- [`feat_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/feat_heatmap.md)
  : Create a heatmap of features organized by cluster
- [`incoming_signaling_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/incoming_signaling_heatmap.md)
  : Create a cluster incoming signaling heatmap
- [`signaling_heatmap()`](https://FertigLab.github.io/dominoSignal/dev/reference/signaling_heatmap.md)
  : Create a network heatmap

### Network Plots

- [`gene_network()`](https://FertigLab.github.io/dominoSignal/dev/reference/gene_network.md)
  : Create a gene association network
- [`signaling_network()`](https://FertigLab.github.io/dominoSignal/dev/reference/signaling_network.md)
  : Create a cluster to cluster signaling network diagram

### Other Plots

- [`circos_ligand_receptor()`](https://FertigLab.github.io/dominoSignal/dev/reference/circos_ligand_receptor.md)
  : Plot expression of a receptor's ligands by other cell types as a
  chord plot
- [`cor_scatter()`](https://FertigLab.github.io/dominoSignal/dev/reference/cor_scatter.md)
  : Create a correlation plot between TF and receptor
- [`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md)
  : Plot differential linkages among domino results ranked by a
  comparative statistic

## Helper Functions

Convenience functions for working with dominoSignal objects

### Access functions

- [`dom_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_clusters.md)
  : Access clusters
- [`dom_correlations()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_correlations.md)
  : Access correlations
- [`dom_counts()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_counts.md)
  : Access counts
- [`dom_database()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_database.md)
  : Access database
- [`dom_de()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_de.md)
  : Access differential expression
- [`dom_info()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_info.md)
  : Access build information
- [`dom_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_linkages.md)
  : Access linkages
- [`dom_network_items()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_network_items.md)
  : Access all features, receptors, or ligands present in a signaling
  network.
- [`dom_signaling()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_signaling.md)
  : Access signaling
- [`dom_tf_activation()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_tf_activation.md)
  : Access transcription factor activation
- [`dom_zscores()`](https://FertigLab.github.io/dominoSignal/dev/reference/dom_zscores.md)
  : Access z-scores

### Misc. utility functions

- [`add_rl_column()`](https://FertigLab.github.io/dominoSignal/dev/reference/add_rl_column.md)
  : Adds a column to the RL signaling data frame.
- [`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
  : Count occurrences of linkages across multiple domino results from a
  linkage summary
- [`mean_ligand_expression()`](https://FertigLab.github.io/dominoSignal/dev/reference/mean_ligand_expression.md)
  : Calculate mean ligand expression as a data frame for plotting in
  circos plot
- [`print(`*`<domino>`*`)`](https://FertigLab.github.io/dominoSignal/dev/reference/print-domino-method.md)
  : Print domino object
- [`print(`*`<linkage_summary>`*`)`](https://FertigLab.github.io/dominoSignal/dev/reference/print-linkage_summary-method.md)
  : Print linkage summary object
- [`rename_clusters()`](https://FertigLab.github.io/dominoSignal/dev/reference/rename_clusters.md)
  : Renames clusters in a domino object
- [`show(`*`<domino>`*`)`](https://FertigLab.github.io/dominoSignal/dev/reference/show-domino-method.md)
  : Show domino object information
- [`show(`*`<linkage_summary>`*`)`](https://FertigLab.github.io/dominoSignal/dev/reference/show-linkage_summary-method.md)
  : Show linkage_summary object information

## Example Data

Example data used for testing and demonstration purposes

- [`CellPhoneDB`](https://FertigLab.github.io/dominoSignal/dev/reference/CellPhoneDB.md)
  : CellPhoneDB subset
- [`SCENIC`](https://FertigLab.github.io/dominoSignal/dev/reference/SCENIC.md)
  : SCENIC AUC subset
- [`PBMC`](https://FertigLab.github.io/dominoSignal/dev/reference/PBMC.md)
  : PBMC RNAseq data subset
- [`mock_linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/mock_linkage_summary.md)
  : Create a mock linkage summary object
