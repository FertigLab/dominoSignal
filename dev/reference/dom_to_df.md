# Turn domino object signaling information into a data frame

Constructs a data frame of ligand-receptor-TF signaling triplets from a
built
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object. Expression values are averaged across cells within each cluster
(ligands from sending clusters, receptors and TFs from receiving
clusters).

## Usage

``` r
dom_to_df(dom, send_clusters = NULL, rec_clusters = NULL, exp_type = "counts")
```

## Arguments

- dom:

  A built
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object (as output by
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

- send_clusters:

  Character/factor vector of cluster names for ligand signals. If NULL
  (default), uses all clusters in dom.

- rec_clusters:

  Character/factor vector of cluster names for receptor/TF signals. If
  NULL (default), uses all clusters in dom.

- exp_type:

  Either "counts" or "z_scores" for expression type

## Value

Data frame with columns: ligand, receptor, transcription_factor,
ligand_exp, rec_exp, tf_auc, sending_cl, receiving_cl. Each row is a
ligand-receptor-TF triplet from a sender-to-receiver cluster pair with
corresponding mean expression values. Empty data frame returned if no
valid interactions found.

## Examples

``` r
data("DominoObjects")
dom <- DominoObjects$built_dom_tiny
df <- dom_to_df(dom, exp_type = "z_scores")
head(df)
#>                  ligand receptor transcription_factor  ligand_exp   rec_exp
#> 1                 CCL20     CCR6                  DBP -0.05812398 0.2108858
#> 2                 CCL20     CCR6                  DBP -0.05812398 0.2108858
#> 3                 CCL20     CCR6                  DBP  0.20959562 0.2108858
#> 4 integrin_a6b4_complex     NRG1                FOSL1  0.11496606 0.4702549
#> 5 integrin_a6b4_complex     NRG1                FOSL1 -0.06554457 0.4702549
#> 6 integrin_a6b4_complex     NRG1                FOSL1 -0.07161490 0.4702549
#>       tf_auc    sending_cl  receiving_cl
#> 1 0.04669679        B_cell        B_cell
#> 2 0.04669679 CD14_monocyte        B_cell
#> 3 0.04669679    CD8_T_cell        B_cell
#> 4 0.06595993        B_cell CD14_monocyte
#> 5 0.06595993 CD14_monocyte CD14_monocyte
#> 6 0.06595993    CD8_T_cell CD14_monocyte
```
