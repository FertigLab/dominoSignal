# Summarize linkages from multiple domino objects

Creates a
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/reference/linkage_summary-class.md)
object storing the linkages learned in different domino objects as
nested lists to facilitate comparisons of networks learned by domino
across subject covariates.

## Usage

``` r
summarize_linkages(domino_results, subject_meta, subject_names = NULL)
```

## Arguments

- domino_results:

  list of domino result with one domino object per subject. Names from
  the list must match subject_names

- subject_meta:

  data frame that includes the subject features by which the objects
  could be grouped. The first column should must be subject names

- subject_names:

  vector of subject names in domino_results. If NULL, defaults to first
  column of subject_meta.

## Value

A linkage summary class object consisting of nested lists of the active
transcription factors, active receptors, and incoming ligands for each
cluster across multiple domino results

## Examples

``` r
example(build_domino, echo = FALSE)

#create alternative clustering by shuffling cluster assignments
clusters_tiny_alt <- setNames(
  PBMC$clusters_tiny[c(121:240, 1:120, 241:360)], 
  names(PBMC$clusters_tiny)
)
clusters_tiny_alt <- as.factor(clusters_tiny_alt)

#build an alternative domino object
pbmc_dom_tiny_alt <- create_domino(
  rl_map = rl_map_tiny,
  features = SCENIC$auc_tiny,
  counts = PBMC$RNA_count_tiny,
  z_scores = PBMC$RNA_zscore_tiny,
  clusters = clusters_tiny_alt,
  tf_targets = regulon_list_tiny,
  use_clusters = TRUE,
  use_complexes = TRUE,
  remove_rec_dropout = FALSE
)
#> Reading in and processing signaling database
#> Database provided from source: CellPhoneDB
#> Getting z_scores, clusters, and counts
#> Calculating feature enrichment by cluster
#> 1 of 3
#> 2 of 3
#> 3 of 3
#> Calculating correlations
#> 1 of 6
#> 2 of 6
#> 3 of 6
#> 4 of 6
#> 5 of 6
#> 6 of 6

pbmc_dom_built_tiny_alt <- build_domino(
  dom = pbmc_dom_tiny_alt,
  min_tf_pval = .05,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .1,
  min_rec_percentage = 0.01
)

#create a list of domino objects
dom_ls <- list(
 dom1 = pbmc_dom_built_tiny,
 dom2 = pbmc_dom_built_tiny_alt
)

#compare the linkages across the two domino objects
meta_df <- data.frame("ID" = c("dom1", "dom2"), "group" = c("A", "B"))
summarize_linkages(
 domino_results = dom_ls, subject_meta = meta_df,
 subject_names = meta_df$ID
)
#> An object of class "linkage_summary"
#> Slot "subject_names":
#> [1] dom1 dom2
#> Levels: dom1 dom2
#> 
#> Slot "subject_meta":
#>     ID group
#> 1 dom1     A
#> 2 dom2     B
#> 
#> Slot "subject_linkages":
#> $dom1
#> $dom1$B_cell
#> $dom1$B_cell$tfs
#> character(0)
#> 
#> $dom1$B_cell$rec
#> NULL
#> 
#> $dom1$B_cell$incoming_lig
#> NULL
#> 
#> $dom1$B_cell$tfs_rec
#> character(0)
#> 
#> $dom1$B_cell$rec_lig
#> character(0)
#> 
#> 
#> $dom1$CD14_monocyte
#> $dom1$CD14_monocyte$tfs
#> [1] "ZNF324" "CREM"   "FOSL1" 
#> 
#> $dom1$CD14_monocyte$rec
#> [1] "CXCR3"        "IL7_receptor" "TGFBR3"       "NRG1"        
#> 
#> $dom1$CD14_monocyte$incoming_lig
#> [1] "CCL20"                 "IL7"                   "TGFB3"                
#> [4] "integrin_a6b4_complex"
#> 
#> $dom1$CD14_monocyte$tfs_rec
#> [1] "ZNF324 <- CXCR3"      "CREM <- IL7_receptor" "CREM <- TGFBR3"      
#> [4] "FOSL1 <- NRG1"       
#> 
#> $dom1$CD14_monocyte$rec_lig
#> [1] "CXCR3 <- CCL20"                "IL7_receptor <- IL7"          
#> [3] "TGFBR3 <- TGFB3"               "NRG1 <- integrin_a6b4_complex"
#> 
#> 
#> $dom1$CD8_T_cell
#> $dom1$CD8_T_cell$tfs
#> [1] "FLI1"
#> 
#> $dom1$CD8_T_cell$rec
#> [1] "CXCR3"        "IL7_receptor"
#> 
#> $dom1$CD8_T_cell$incoming_lig
#> [1] "CCL20" "IL7"  
#> 
#> $dom1$CD8_T_cell$tfs_rec
#> [1] "FLI1 <- CXCR3"        "FLI1 <- IL7_receptor"
#> 
#> $dom1$CD8_T_cell$rec_lig
#> [1] "CXCR3 <- CCL20"      "IL7_receptor <- IL7"
#> 
#> 
#> 
#> $dom2
#> $dom2$B_cell
#> $dom2$B_cell$tfs
#> character(0)
#> 
#> $dom2$B_cell$rec
#> NULL
#> 
#> $dom2$B_cell$incoming_lig
#> NULL
#> 
#> $dom2$B_cell$tfs_rec
#> character(0)
#> 
#> $dom2$B_cell$rec_lig
#> character(0)
#> 
#> 
#> $dom2$CD14_monocyte
#> $dom2$CD14_monocyte$tfs
#> [1] "FLI1"
#> 
#> $dom2$CD14_monocyte$rec
#> [1] "CXCR3"        "IL7_receptor"
#> 
#> $dom2$CD14_monocyte$incoming_lig
#> [1] "CCL20" "IL7"  
#> 
#> $dom2$CD14_monocyte$tfs_rec
#> [1] "FLI1 <- CXCR3"        "FLI1 <- IL7_receptor"
#> 
#> $dom2$CD14_monocyte$rec_lig
#> [1] "CXCR3 <- CCL20"      "IL7_receptor <- IL7"
#> 
#> 
#> $dom2$CD8_T_cell
#> $dom2$CD8_T_cell$tfs
#> [1] "ZNF324" "CREM"   "FOSL1" 
#> 
#> $dom2$CD8_T_cell$rec
#> [1] "CXCR3"        "IL7_receptor" "TGFBR3"       "NRG1"        
#> 
#> $dom2$CD8_T_cell$incoming_lig
#> [1] "CCL20"                 "IL7"                   "TGFB3"                
#> [4] "integrin_a6b4_complex"
#> 
#> $dom2$CD8_T_cell$tfs_rec
#> [1] "ZNF324 <- CXCR3"      "CREM <- IL7_receptor" "CREM <- TGFBR3"      
#> [4] "FOSL1 <- NRG1"       
#> 
#> $dom2$CD8_T_cell$rec_lig
#> [1] "CXCR3 <- CCL20"                "IL7_receptor <- IL7"          
#> [3] "TGFBR3 <- TGFB3"               "NRG1 <- integrin_a6b4_complex"
#> 
#> 
#> 
#> 
```
