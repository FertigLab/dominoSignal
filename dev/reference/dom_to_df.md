# Turn domino object signaling information into a data frame

Constructs a data frame of ligand-receptor-TF signaling triplets from a
built domino object. Expression values are averaged across cells within
each cluster (ligands from sending clusters, receptors and TFs from
receiving clusters).

## Usage

``` r
dom_to_df(
  dom,
  send_clusters = NULL,
  rec_clusters = NULL,
  exp_type = c("counts", "z_scores")
)
```

## Arguments

- dom:

  A built domino object (as output by
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md))

- send_clusters:

  Character/factor vector of cluster names for ligand signals. If NULL
  (default), uses all clusters in dom.

- rec_clusters:

  Character/factor vector of cluster names for receptor/TF signals. If
  NULL (default), uses all clusters in dom.

- exp_type:

  Character of length 1: either "counts" or "z_scores" for expression
  type

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
#>   ligand     receptor transcription_factor  ligand_exp     rec_exp     tf_auc
#> 1  CCL20        CXCR3               ZNF324 -0.05812398 -0.05510291 0.03354529
#> 2  CCL20        CXCR3               ZNF324  0.20959562 -0.05510291 0.03354529
#> 3  CCL20        CXCR3               ZNF324 -0.05812398 -0.05510291 0.03354529
#> 4    IL7 IL7_receptor                 CREM  0.16609099 -0.24893817 0.04772885
#> 5    IL7 IL7_receptor                 CREM -0.08080662 -0.24893817 0.04772885
#> 6    IL7 IL7_receptor                 CREM  0.08364176 -0.24893817 0.04772885
#>      sending_cl  receiving_cl
#> 1        B_cell CD14_monocyte
#> 2 CD14_monocyte CD14_monocyte
#> 3    CD8_T_cell CD14_monocyte
#> 4        B_cell CD14_monocyte
#> 5 CD14_monocyte CD14_monocyte
#> 6    CD8_T_cell CD14_monocyte
```
