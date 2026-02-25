# Plotting Functions and Options

In addition to the numerous plotting options available in the original
[Domino](https://github.com/Elisseeff-Lab/domino), dominoSignal has
added more functionality and new methods to improve visualization and
interpretation of analysis results. Here, we will go over the different
plotting functions available, as well as different options that can be
utilized, and options for customization.

## Setup and Data Load

``` r

set.seed(42)
library(dominoSignal)
library(patchwork)
```

In this tutorial, we will use the domino object we built on the [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page. If you have not yet built a domino object, you can do so by
following the instructions on the [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page.

Instructions to load data

``` r

# BiocFileCache helps with managing files across sessions
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
data_url <- "https://zenodo.org/records/10951634/files/pbmc_domino_built.rds"
tmp_path <- BiocFileCache::bfcrpath(bfc, data_url)
dom <- readRDS(tmp_path)
```

## Heatmaps

### Correlations between receptors and transcription factors

[`cor_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/cor_heatmap.md)
can be used to show the correlations calculated between receptors and
transcription factors.

``` r

cor_heatmap(dom, title = "PBMC R-TF Correlations", column_names_gp = grid::gpar(fontsize = 8))
```

![Heatmap of correlation values between receptors and transcription
factors across the
set.](plotting_vignette_files/figure-html/corheatmap-1.png)

In addition to displaying the scores for the correlations, this function
can also be used to identify correlations above a certain value (using
the `bool` and `bool_thresh` arguments) or to identify the combinations
of receptors and transcription factors (TFs) that are connected (with
argument `mark_connections`).

``` r

cor_heatmap(dom, bool = TRUE, bool_thresh = 0.25)
cor_heatmap(dom, bool = FALSE, mark_connections = TRUE)
```

![Heatmap of correlations shown as true or false with threshold of 0.25
and heatmap of correlations with xs shown over receptors and
transcription factors that are
linked.](plotting_vignette_files/figure-html/corheatmap-options-1.png)![Heatmap
of correlations shown as true or false with threshold of 0.25 and
heatmap of correlations with xs shown over receptors and transcription
factors that are
linked.](plotting_vignette_files/figure-html/corheatmap-options-2.png)

If only a subset of receptors or transcription factors are of interest,
a vector of either (or both) can be passed to the function.

``` r

receptors <- c("CSF1R", "CSF3R", "CCR7", "FCER2")
tfs <- c("PAX5", "JUNB", "FOXJ3", "FOSB")
cor_heatmap(dom, feats = tfs, recs = receptors)
```

![A subsetted heatmap showing the correlation between specific
transcription factors and receptors of
interest.](plotting_vignette_files/figure-html/corheatmap-subset-1.png)

The heatmap functions in dominoSignal are based on
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
and will also accept additional arguments meant for that function. For
example, while an argument for clustering rows or columns of the heatmap
is not explicitly stated, they can still be passed to
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
through
[`cor_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/cor_heatmap.md).

``` r

cor_heatmap(dom, cluster_rows = FALSE, cluster_columns = FALSE, column_title = "Heatmap Without Clustering",
    column_names_gp = grid::gpar(fontsize = 4))
```

![A correlation heatmap that has not been clustered and does not display
a dendrogram as a
result.](plotting_vignette_files/figure-html/corheatmap-ComplexHeatmap-args-1.png)

### Heatmap of Transcription Factor Activation Scores

[`feat_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/feat_heatmap.md)
is used to show the transcription factor activation for features in the
signaling network.

``` r

feat_heatmap(dom, use_raster = FALSE, row_names_gp = grid::gpar(fontsize = 4))
```

![A heatmap showing transcription factor activation levels by cells
which are arranged by
cluster.](plotting_vignette_files/figure-html/featheatmap-1.png)

It functions similarly to
[`cor_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/cor_heatmap.md),
with arguments to select a specific vector of features, to use a boolean
view with a threshold, and to pass other arguments to
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).
Specific to this function though are arguments for setting the range of
values to be visualized and one to choose to normalize the scores to the
max value.

``` r

feat_heatmap(dom, min_thresh = 0.1, max_thresh = 0.6, norm = TRUE, bool = FALSE,
    use_raster = FALSE)
feat_heatmap(dom, bool = TRUE, use_raster = FALSE)
```

![A heatmap showing activation of transcription factors with a minimum
activity of 0.1 and maximum activity of 0.6 with values normalized to
the maximum score, and a binary heatmap showing transcription factor
activation with the default threshold of
0.2.](plotting_vignette_files/figure-html/featheatmap-options-1.png)![A
heatmap showing activation of transcription factors with a minimum
activity of 0.1 and maximum activity of 0.6 with values normalized to
the maximum score, and a binary heatmap showing transcription factor
activation with the default threshold of
0.2.](plotting_vignette_files/figure-html/featheatmap-options-2.png)

### Heatmap of Incoming Signaling for a Cluster

[`incoming_signaling_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/incoming_signaling_heatmap.md)
can be used to visualize the cluster average expression of the ligands
capable of activating the TFs enriched in the cluster. For example, to
view the incoming signaling of the CD8 T cells:

``` r

incoming_signaling_heatmap(dom, "CD8_T_cell")
```

![Heatmap of ligand expression by sending cluster targetting CD8 T
cells.](plotting_vignette_files/figure-html/incoming-1.png)

We can also select for specific clusters of interest that are signaling
to the CD8 T cells. If we are only interested in viewing the monocyte
signaling:

``` r

incoming_signaling_heatmap(dom, "CD8_T_cell", clusts = c("CD14_monocyte", "CD16_monocyte"))
```

![Heatmap of ligand expression subset to sending clusters CD14 monocytes
and CD16 monocytes targetting CD8 T
cells.](plotting_vignette_files/figure-html/incoming-subset-1.png)

As with our other heatmap functions, options are available for a minimum
threshold, maximum threshold, whether to scale the values after
thresholding, whether to normalize the matrix, and the ability to pass
further arguments to
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

### Heatmap of Signaling Between Clusters

[`signaling_heatmap()`](https://FertigLab.github.io/dominoSignal/reference/signaling_heatmap.md)
makes a heatmap showing the signaling strength of ligands from each
cluster to receptors of each cluster based on averaged expression.

``` r

signaling_heatmap(dom)
```

![Heatmap of collective signaling between receiving clusters and sending
clusters.](plotting_vignette_files/figure-html/signaling-1.png)

As with other functions, specific clusters can be selected, thresholds
can be set, and normalization methods can be selected as well.

``` r

signaling_heatmap(dom, scale = "sqrt")
signaling_heatmap(dom, normalize = "rec_norm")
```

![A heatmap of collective signaling with a square root transformation
and a heatmap of collective signaling normalized to the maximum value of
received
signaling.](plotting_vignette_files/figure-html/signaling-norm-1.png)![A
heatmap of collective signaling with a square root transformation and a
heatmap of collective signaling normalized to the maximum value of
received
signaling.](plotting_vignette_files/figure-html/signaling-norm-2.png)

## Network Plots

### Network showing L - R - TF signaling between clusters

[`gene_network()`](https://FertigLab.github.io/dominoSignal/reference/gene_network.md)
makes a network plot to display signaling in selected clusters such that
ligands, receptors and features associated with those clusters are
displayed as nodes with edges as linkages. To look at signaling to the
CD16 Monocytes from the CD14 Monocytes:

``` r

# \ linkages from the CD14 monocyte cluster to the CD16 monocyte cluster.
gene_network(dom, clust = "CD16_monocyte", OutgoingSignalingClust = "CD14_monocyte")
```

![Network plot with columns for ligand, receptor, and transcription
factor showing](plotting_vignette_files/figure-html/genenetwork-1.png)

Options to modify this plot include adjusting scaling for the ligands
and different layouts (some of which are more legible than others).

``` r

gene_network(dom, clust = "CD16_monocyte", OutgoingSignalingClust = "CD14_monocyte",
    lig_scale = 25, layout = "sphere")
```

![Network plot for CD14 monocytes to CD16 monocyte with ligand values
scaled up and a spherical network
layout.](plotting_vignette_files/figure-html/genenetwork-options-1.png)

Additionally, colors can be given for select genes (for example, to
highlight a specific signaling path).

``` r

gene_network(dom, clust = "CD16_monocyte", OutgoingSignalingClust = "CD14_monocyte",
    cols = c(CD1D = "violet", LILRB2 = "violet", FOSB = "violet"), lig_scale = 10)
```

![Network plot from CD14 monocytes to CD16 monocytes in three column
layout where ligand CD1D, receptor LILRB2, and transcription factor FOSB
have been highlighted in violet and ligand nodes are scaled
up.](plotting_vignette_files/figure-html/genenetwork-cols-1.png)

### Network Showing Interaction Strength Across Data

[`signaling_network()`](https://FertigLab.github.io/dominoSignal/reference/signaling_network.md)
can be used to create a network plot such that nodes are clusters and
the edges indicate signaling from one cluster to another.

``` r

signaling_network(dom)
```

![Network plot with circular layout with nodes corresponding to clusters
and edges showing connections between clusters. The edges from the
platelet population dominate the
plot.](plotting_vignette_files/figure-html/signalingnet-1.png)

In addition to changes such as scaling, thresholds, or layouts, this
plot can be modified by selecting incoming or outgoing clusters (or
both!). An example to view signaling from the CD14 Monocytes to other
clusters:

``` r

signaling_network(dom, showOutgoingSignalingClusts = "CD14_monocyte", scale = "none",
    norm = "none", layout = "fr", scale_by = "none", edge_weight = 2, vert_scale = 10)
```

![Network plot for signaling from the CD14 monocytes which are plotted
in the center with edges connecting to the other clusters which are
positioned around the central
node.](plotting_vignette_files/figure-html/signalingnet-clusts-1.png)

## Other Types of Plots

### Chord Diagrams Connecting Ligands and Receptors

A new addition to dominoSignal,
[`circos_ligand_receptor()`](https://FertigLab.github.io/dominoSignal/reference/circos_ligand_receptor.md)
creates a chord plot showing ligands that can activate a specific
receptor, displaying mean cluster expression of the ligand with the
width of the chord.

``` r

circos_ligand_receptor(dom, receptor = "CD74")
```

![Circos plot with outer arcs corresponding to sending clusters and
inner arcs corresponding to ligands connecting to receptor
CD74.](plotting_vignette_files/figure-html/circos-1.png)

This function can be given cluster colors to match other plots you may
make with your data. In addition, the plot can be adjusted by changing
the threshold of ligand expression required for a linkage to be
visualized or selecting clusters of interest.

``` r

cols <- c("red", "orange", "green", "blue", "pink", "purple", "slategrey", "firebrick",
    "hotpink")
names(cols) <- dom_clusters(dom, labels = FALSE)
circos_ligand_receptor(dom, receptor = "CD74", cell_colors = cols)
```

![Circos plot identical to the one above, except that the outer cluster
arcs have been assigned specific
colors.](plotting_vignette_files/figure-html/circos-opt-1.png)

### Scatter Plot to Visualize Correlation

[`cor_scatter()`](https://FertigLab.github.io/dominoSignal/reference/cor_scatter.md)
can be used to plot each cell based on activation of the selected TF and
expression of the receptor. This produces a scatter plot as well as a
line of best fit to look at receptor - TF correlation.

``` r

cor_scatter(dom, "FOSB", "CD74")
```

![Scatter plot of receptor expression by transcription factor activation
with a line of best fit displayed over the
points.](plotting_vignette_files/figure-html/corscatter-1.png)

Do keep in mind that there is an argument for `remove_rec_dropout` that
should match the parameter that was used when the domino object was
built. In this case, we did not use that build parameter, so we will
leave the value in this argument as its default value of FALSE.

## Continued Development

Since dominoSignal is a package still being developed, there are new
functions and features that will be implemented in future versions. In
the meantime, to view an example analysis, see our [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page, or see our [domino object
structure](https://fertiglab.github.io/dominoSignal/articles/domino_object_vignette)
page to get familiar with the object structure. Additionally, if you
find any bugs, have further questions, or want to share an idea, please
let us know [here](https://github.com/FertigLab/dominoSignal/issues).

Vignette Build Information

Date last built and session information:

``` r

Sys.Date()
#> [1] "2026-02-25"
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.2
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] patchwork_1.3.2    dominoSignal_1.0.5
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.2.3             httr2_1.2.1           formatR_1.14         
#>   [4] biomaRt_2.66.0        rlang_1.1.6           magrittr_2.0.4       
#>   [7] clue_0.3-66           GetoptLong_1.1.0      matrixStats_1.5.0    
#>  [10] compiler_4.5.2        RSQLite_2.4.5         mgcv_1.9-4           
#>  [13] png_0.1-8             systemfonts_1.3.1     vctrs_0.6.5          
#>  [16] stringr_1.6.0         pkgconfig_2.0.3       shape_1.4.6.1        
#>  [19] crayon_1.5.3          fastmap_1.2.0         backports_1.5.0      
#>  [22] dbplyr_2.5.1          XVector_0.50.0        labeling_0.4.3       
#>  [25] rmarkdown_2.30        ragg_1.5.0            purrr_1.2.0          
#>  [28] bit_4.6.0             xfun_0.54             cachem_1.1.0         
#>  [31] jsonlite_2.0.0        progress_1.2.3        blob_1.2.4           
#>  [34] broom_1.0.11          parallel_4.5.2        prettyunits_1.2.0    
#>  [37] cluster_2.1.8.1       R6_2.6.1              bslib_0.9.0          
#>  [40] stringi_1.8.7         RColorBrewer_1.1-3    car_3.1-3            
#>  [43] jquerylib_0.1.4       Rcpp_1.1.0            Seqinfo_1.0.0        
#>  [46] iterators_1.0.14      knitr_1.50            IRanges_2.44.0       
#>  [49] splines_4.5.2         Matrix_1.7-4          igraph_2.2.1         
#>  [52] tidyselect_1.2.1      dichromat_2.0-0.1     abind_1.4-8          
#>  [55] yaml_2.3.11           doParallel_1.0.17     codetools_0.2-20     
#>  [58] curl_7.0.0            lattice_0.22-7        tibble_3.3.0         
#>  [61] plyr_1.8.9            withr_3.0.2           Biobase_2.70.0       
#>  [64] KEGGREST_1.50.0       S7_0.2.1              evaluate_1.0.5       
#>  [67] desc_1.4.3            BiocFileCache_3.0.0   circlize_0.4.16      
#>  [70] Biostrings_2.78.0     pillar_1.11.1         ggpubr_0.6.2         
#>  [73] filelock_1.0.3        carData_3.0-5         foreach_1.5.2        
#>  [76] stats4_4.5.2          generics_0.1.4        S4Vectors_0.48.0     
#>  [79] hms_1.1.4             ggplot2_4.0.1         scales_1.4.0         
#>  [82] glue_1.8.0            tools_4.5.2           ggsignif_0.6.4       
#>  [85] fs_1.6.6              Cairo_1.7-0           grid_4.5.2           
#>  [88] tidyr_1.3.1           AnnotationDbi_1.72.0  colorspace_2.1-2     
#>  [91] nlme_3.1-168          Formula_1.2-5         cli_3.6.5            
#>  [94] rappdirs_0.3.3        textshaping_1.0.4     ComplexHeatmap_2.26.0
#>  [97] dplyr_1.1.4           gtable_0.3.6          rstatix_0.7.3        
#> [100] sass_0.4.10           digest_0.6.39         BiocGenerics_0.56.0  
#> [103] rjson_0.2.23          htmlwidgets_1.6.4     farver_2.1.2         
#> [106] memoise_2.0.1         htmltools_0.5.9       pkgdown_2.2.0        
#> [109] lifecycle_1.0.4       httr_1.4.7            GlobalOptions_0.1.3  
#> [112] bit64_4.6.0-1
```
