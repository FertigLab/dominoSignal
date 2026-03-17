# Interacting with domino Objects

In this tutorial, we will go into more detail on the structure of a
domino object and the ways in which to access the data stored within. We
will be using the domino object we built on the [Getting
Started](https://fertiglab.github.io/dominoSignal) page. If you have not
yet built a domino object, you can do so by following the instructions
on the [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page.

Instructions to load data

``` r

set.seed(42)
library(dominoSignal)

# BiocFileCache helps with file management across sessions
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
data_url <- "https://zenodo.org/records/10951634/files/pbmc_domino_built.rds"
tmp_path <- BiocFileCache::bfcrpath(bfc, data_url)
dom <- readRDS(tmp_path)
```

## Object contents

There is a great deal of information stored with the [domino object
class](https://fertiglab.github.io/dominoSignal/reference/domino-class.html).
The domino object is an S4 class object that contains a variety of
information about the data set used to build the object, the calculated
values, and the linkages between receptors, ligands, and transcription
factors. The object is structured as follows (with some examples of the
information stored within each slot:

- Input Data

  - Information about the database used to construct the rl_map

  - Inputted counts matrix

  - Inputted z-scored counts matrix

  - Inputted cluster labels

  - Inputted transcription factor activation scores

- Calculated values

  - Differential expression p-values of transcription factors in each
    cluster

  - Correlation values between ligands and receptors

  - Median correlation between components of receptor complexes

- Linkages

  - Complexes show the component genes of any complexes in the rl map

  - Receptor - ligand linkages as determined from the rl map

  - Transcription factor - target linkages as determined from the SCENIC
    analysis (or other regulon inference method)

  - Transcription factors that are differentially expressed in each
    cluster

  - Transcription factors that are correlated with receptors

  - Transcription factors that are correlated with receptors in each
    cluster

  - Receptors which are active in each cluster

  - Ligands that may activate a receptor in a given cluster (so-called
    incoming ligands; these may include ligands from outside the data
    set)

- Signaling matrices

  - For each cluster, incoming ligands and the clusters within the data
    set that they are coming from

  - A summary of signaling between all clusters

- Miscellaneous Information

  - Build information, which includes the parameters used to build the
    object in the
    [`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md)
    functions

  - The pared down receptor ligand map information used in building the
    object

  - The percent expression of receptors within each cluster

For commonly accessed information (the number of cells, clusters, and
some build information), the show and print methods for domino objects
can be used.

``` r

dom
#> A domino object of 2607 cells
#> Built with signaling between 9 clusters
```

``` r

print(dom)
```

## Access functions

To facilitate access to the information stored in the domino object, we
have provided a collection of functions to retrieve specific items.
These functions begin with “dom\_” and can be listed using
[`ls()`](https://rdrr.io/r/base/ls.html).

``` r

ls("package:dominoSignal", pattern = "^dom_")
#>  [1] "dom_clusters"      "dom_correlations"  "dom_counts"       
#>  [4] "dom_database"      "dom_de"            "dom_info"         
#>  [7] "dom_linkages"      "dom_network_items" "dom_signaling"    
#> [10] "dom_tf_activation" "dom_zscores"
```

### Input data

When creating a domino object with the
[`create_domino()`](https://FertigLab.github.io/dominoSignal/reference/create_domino.md)
function, several inputs are required which are then stored in the
domino object itself. These include cluster labels, the counts matrix,
z-scored counts matrix, transcription factor activation scores, and the
R-L database used in
[`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/reference/create_rl_map_cellphonedb.md).

For example, to access the cluster names in the domino object:

``` r

dom_clusters(dom)
#> [1] "B_cell"            "CD14_monocyte"     "CD16_monocyte"    
#> [4] "CD8_T_cell"        "dendritic_cell"    "memory_CD4_T_cell"
#> [7] "naive_CD4_T_cell"  "NK_cell"           "Platelet"
```

Setting an argument `labels = TRUE` will return the vector of cluster
labels for each cell rather than the unique cluster names.

To access the counts:

``` r

count_matrix <- dom_counts(dom)
knitr::kable(count_matrix[1:5, 1:5])
```

|  | AAACATACAACCAC-1 | AAACATTGAGCTAC-1 | AAACATTGATCAGC-1 | AAACCGTGCTTCCG-1 | AAACCGTGTATGCG-1 |
|:---|---:|---:|---:|---:|---:|
| AL627309.1 | 0 | 0 | 0 | 0 | 0 |
| AP006222.2 | 0 | 0 | 0 | 0 | 0 |
| RP11-206L10.2 | 0 | 0 | 0 | 0 | 0 |
| RP11-206L10.9 | 0 | 0 | 0 | 0 | 0 |
| FAM87B | 0 | 0 | 0 | 0 | 0 |

Or z-scored counts:

``` r

z_matrix <- dom_zscores(dom)
knitr::kable(z_matrix[1:5, 1:5])
```

|  | AAACATACAACCAC-1 | AAACATTGAGCTAC-1 | AAACATTGATCAGC-1 | AAACCGTGCTTCCG-1 | AAACCGTGTATGCG-1 |
|:---|---:|---:|---:|---:|---:|
| AL627309.1 | -0.0581122 | -0.0581122 | -0.0581122 | -0.0581122 | -0.0581122 |
| AP006222.2 | -0.0335151 | -0.0335151 | -0.0335151 | -0.0335151 | -0.0335151 |
| RP11-206L10.2 | -0.0399375 | -0.0399375 | -0.0399375 | -0.0399375 | -0.0399375 |
| RP11-206L10.9 | -0.0337574 | -0.0337574 | -0.0337574 | -0.0337574 | -0.0337574 |
| FAM87B | -0.0274227 | -0.0274227 | -0.0274227 | -0.0274227 | -0.0274227 |

The transcription factor activation scores can be similarly accessed:

``` r

activation_matrix <- dom_tf_activation(dom)
knitr::kable(activation_matrix[1:5, 1:5])
```

|  | AAACATACAACCAC-1 | AAACATTGAGCTAC-1 | AAACATTGATCAGC-1 | AAACCGTGCTTCCG-1 | AAACCGTGTATGCG-1 |
|:---|---:|---:|---:|---:|---:|
| ARNTL | 0.0446386 | 0.0413253 | 0.0441566 | 0.0437952 | 0.1087952 |
| ATF3 | 0.0582029 | 0.0870623 | 0.0509289 | 0.1172953 | 0.0638104 |
| ATF4 | 0.0295733 | 0.0437047 | 0.0542889 | 0.0427879 | 0.0436855 |
| ATF6 | 0.0000000 | 0.0211198 | 0.1151665 | 0.0850461 | 0.0460666 |
| BCL3 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 |

Information about the database referenced for ligand - receptor pairs
and composition of protein complexes can be extracted from the
[`dom_database()`](https://FertigLab.github.io/dominoSignal/reference/dom_database.md)
function. By default, the function returns the name(s) of the
database(s) used:

``` r

dom_database(dom)
#> [1] "CellPhoneDB_v4.0"
```

If you would like to view the entire ligand - receptor map, set
`name_only = FALSE`:

``` r

db_matrix <- dom_database(dom, name_only = FALSE)
knitr::kable(db_matrix[1:5, 1:5])
```

|  | int_pair | name_A | uniprot_A | gene_A | type_A |
|:---|:---|:---|:---|:---|:---|
| 4 | RAreceptor_RXRG & atRetinoicAcid_byALDH1A3 | RAreceptor_RXRG | P48443,P29373 | RXRG,CRABP2 | R |
| 5 | RAreceptor_RXRG & atRetinoicAcid_byALDH1A2 | RAreceptor_RXRG | P48443,P29373 | RXRG,CRABP2 | R |
| 6 | RAreceptor_RXRG & atRetinoicAcid_byALDH1A1 | RAreceptor_RXRG | P48443,P29373 | RXRG,CRABP2 | R |
| 7 | RAreceptor_RXRB & atRetinoicAcid_byALDH1A3 | RAreceptor_RXRB | P28702,P29373 | RXRB,CRABP2 | R |
| 8 | RAreceptor_RXRB & atRetinoicAcid_byALDH1A2 | RAreceptor_RXRB | P28702,P29373 | RXRB,CRABP2 | R |

### Calculations

Active transcription factors in each cluster are determined by
conducting Wilcoxon rank sum tests for each transcription factor where
the transcription factor activity scores amongst all cells in a cluster
are tested against the activity scores of all cells outside of the
cluster. The p-values for the one-sided test for greater activity within
the cluster compared to other cells can be accessed with the
[`dom_de()`](https://FertigLab.github.io/dominoSignal/reference/dom_de.md)
function.

``` r

de_matrix <- dom_de(dom)
knitr::kable(de_matrix[1:5, 1:5])
```

|       |    B_cell | CD14_monocyte | CD16_monocyte | CD8_T_cell | dendritic_cell |
|:------|----------:|--------------:|--------------:|-----------:|---------------:|
| ARNTL | 0.0000019 |             1 |     0.9947780 |  0.9865518 |      0.8359381 |
| ATF3  | 0.0006972 |             0 |     0.0000000 |  1.0000000 |      0.0000000 |
| ATF4  | 0.8399818 |             1 |     1.0000000 |  0.9631934 |      0.1748338 |
| ATF6  | 0.5988625 |             1 |     0.9893142 |  0.0365779 |      0.9549489 |
| BCL3  | 0.9299696 |             0 |     0.6007326 |  0.1482247 |      0.0757460 |

Linkage between receptors and transcription factors is assessed by
Spearman correlation between transcription factor activity scores and
scaled expression of receptor-encoding genes across all cells in the
data set. Spearman coefficients can be accessed with the
[`dom_correlations()`](https://FertigLab.github.io/dominoSignal/reference/dom_correlations.md)
function. Setting `type` to “complex” will return the median correlation
between components of receptor complexes; the default (“rl”) will return
receptor - ligand correlations.

``` r

cor_matrix <- dom_correlations(dom)
knitr::kable(cor_matrix[1:5, 1:5])
```

|          |      ARNTL |       ATF3 |       ATF4 |       ATF6 |       BCL3 |
|:---------|-----------:|-----------:|-----------:|-----------:|-----------:|
| TNFRSF18 |  0.0057014 | -0.0297120 |  0.0043825 |  0.0063055 | -0.0029035 |
| TNFRSF4  |  0.0373481 | -0.1012705 |  0.0344042 |  0.0508720 |  0.0282500 |
| TNFRSF14 | -0.0356892 |  0.0541435 |  0.0015811 |  0.0093978 |  0.0206700 |
| TNFRSF25 |  0.0310362 | -0.1171502 |  0.0701541 |  0.0458227 |  0.0047186 |
| TNFRSF1B | -0.0933019 |  0.2690328 | -0.0625319 | -0.0457366 |  0.0239629 |

### Linkages

Linkages between ligands, receptors, and transcription factors can be
accessed in several different ways, depending on the specific link and
the scope desired. The
[`dom_linkages()`](https://FertigLab.github.io/dominoSignal/reference/dom_linkages.md)
function has three arguments - the first, like all of our access
functions, is for the domino object. The second, `link_type`, is used to
specify which linkages are desired (options are complexes, receptor -
ligand, tf - target, or tf - receptor). The third argument,
`by_cluster`, determines whether the linkages returned are arranged by
cluster (though this does change the available linkage types to tf -
receptor, receptor, or incoming-ligand). For example, to access the
complexes used across the dataset:

``` r

complex_links <- dom_linkages(dom, link_type = "complexes")
# Look for components of NODAL receptor complex
complex_links$NODAL_receptor
#> NULL
```

To view incoming ligands to each cluster:

``` r

incoming_links <- dom_linkages(dom, link_type = "incoming-ligand", by_cluster = TRUE)
# Check incoming signals to dendritic cells
incoming_links$dendritic_cell
#>  [1] "COPA"                  "MIF"                   "APP"                  
#>  [4] "FAM19A4"               "TAFA4"                 "ANXA1"                
#>  [7] "CD99"                  "integrin_aVb3_complex" "integrin_a4b1_complex"
#> [10] "BMP8B"                 "PLAU"                  "CSF3"                 
#> [13] "CXCL9"                 "HLA-F"                 "CD1D"                 
#> [16] "INS"                   "IL34"                  "CSF1"                 
#> [19] "CSF2"                  "CTLA4"                 "CD28"                 
#> [22] "GRN"                   "TNF"                   "LTA"                  
#> [25] "CXCL12"                "CXCL14"                "CD58"
```

If, for some reason, you find yourself in need of the entire linkage
structure (not recommended), it can be accessed through its slot name;
domino objects are S4 objects.

``` r

all_linkages <- slot(dom, "linkages")
# Names of all sub-structures:
names(all_linkages)
#> [1] "complexes"          "rec_lig"            "tf_targets"        
#> [4] "clust_tf"           "tf_rec"             "clust_tf_rec"      
#> [7] "clust_rec"          "clust_incoming_lig"
```

Alternately, to obtain a simplified list of receptors, ligands, and/or
features in the domino object, use the
[`dom_network_items()`](https://FertigLab.github.io/dominoSignal/reference/dom_network_items.md)
function. To pull all transcription factors associated with the
dendritic cell cluster:

``` r

dc_tfs <- dom_network_items(dom, "dendritic_cell", return = "features")
head(dc_tfs)
#> [1] "ATF3"  "CEBPD" "FOSB"  "STAT6" "KLF4"  "CEBPA"
```

### Signaling Matrices

The averaged z-scored expression of ligands and receptors between
different clusters can be accessed in matrix form.

``` r

signal_matrix <- dom_signaling(dom)
knitr::kable(signal_matrix)
```

|  | L_B_cell | L_CD14_monocyte | L_CD16_monocyte | L_CD8_T_cell | L_dendritic_cell | L_memory_CD4_T_cell | L_naive_CD4_T_cell | L_NK_cell | L_Platelet |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| R_B_cell | 0.0494729 | 0.1976263 | 0.1686996 | 0.0000000 | 0.3040989 | 0.2459618 | 0.1367747 | 0.0000000 | 4.013299 |
| R_CD14_monocyte | 0.8060665 | 2.5378490 | 1.4429140 | 1.1787148 | 3.4961650 | 2.2611486 | 1.2608492 | 1.4485708 | 10.559449 |
| R_CD16_monocyte | 0.8060665 | 2.5378490 | 1.1473216 | 1.1787148 | 3.4961650 | 2.2611486 | 0.8757246 | 1.4055805 | 7.258572 |
| R_CD8_T_cell | 0.3759111 | 0.0000000 | 0.1263349 | 0.2120444 | 0.0000000 | 0.2772703 | 0.0000000 | 0.1604005 | 0.000000 |
| R_dendritic_cell | 0.3176857 | 2.5041485 | 0.7368126 | 1.3152485 | 3.4961650 | 2.3809976 | 0.3061580 | 1.5659811 | 7.258572 |
| R_memory_CD4_T_cell | 0.3759111 | 0.0000000 | 0.1263349 | 0.2120444 | 0.0000000 | 0.3801282 | 0.0000000 | 0.1604005 | 0.000000 |
| R_naive_CD4_T_cell | 0.3759111 | 0.0000000 | 0.1263349 | 0.2120444 | 0.0000000 | 0.3801282 | 0.0000000 | 0.1604005 | 0.000000 |
| R_NK_cell | 0.2930600 | 1.5760509 | 0.6086536 | 0.2556741 | 1.5025648 | 0.7620921 | 0.1772729 | 0.1604005 | 4.013299 |
| R_Platelet | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.000000 |

To view signaling to a specific cluster from the other clusters, set the
`cluster` argument to the cluster name.

``` r

dc_matrix <- dom_signaling(dom, "dendritic_cell")
knitr::kable(dc_matrix)
```

|  | dendritic_cell.L_B_cell | dendritic_cell.L_CD14_monocyte | dendritic_cell.L_CD16_monocyte | dendritic_cell.L_CD8_T_cell | dendritic_cell.L_dendritic_cell | dendritic_cell.L_memory_CD4_T_cell | dendritic_cell.L_naive_CD4_T_cell | dendritic_cell.L_NK_cell | dendritic_cell.L_Platelet |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COPA | 0.0000000 | 0.1352458 | 0.1686996 | 0.0000000 | 0.1706852 | 0.0243003 | 0.0000000 | 0.0000000 | 0.1634935 |
| MIF | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0810828 | 0.1637432 | 0.1367747 | 0.0000000 | 0.0000000 |
| APP | 0.0494729 | 0.0623805 | 0.0000000 | 0.0000000 | 0.0523309 | 0.0579182 | 0.0000000 | 0.0000000 | 3.8498060 |
| ANXA1 | 0.0000000 | 0.1864958 | 0.0000000 | 0.1682595 | 0.6876552 | 0.5065038 | 0.0000000 | 0.3542401 | 0.0000000 |
| CD99 | 0.0000000 | 0.0626658 | 0.0000000 | 0.4718174 | 0.0000000 | 0.1952307 | 0.0000000 | 0.7924524 | 1.9983051 |
| integrin_a4b1_complex | 0.0775377 | 0.0000000 | 0.0274081 | 0.0128595 | 0.0933263 | 0.1408531 | 0.0000000 | 0.0000000 | 1.2469674 |
| BMP8B | 0.1291568 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 |
| PLAU | 0.0000000 | 0.0455798 | 0.0000000 | 0.0711901 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 |
| CXCL9 | 0.0000000 | 0.0663249 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0582573 | 0.0000000 | 0.0000000 | 0.0000000 |
| HLA-F | 0.0000000 | 0.0000000 | 0.0000000 | 0.2918047 | 0.0000000 | 0.0000000 | 0.0000000 | 0.2588880 | 0.0000000 |
| CD1D | 0.0615183 | 0.5670314 | 0.1007509 | 0.0000000 | 1.2126187 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 |
| CSF1 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0067252 | 0.0000000 | 0.1724684 | 0.0000000 | 0.0000000 | 0.0000000 |
| CTLA4 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0369179 | 0.0000000 | 0.2085753 | 0.0000000 | 0.0000000 | 0.0000000 |
| CD28 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.2771802 | 0.1288850 | 0.0000000 | 0.0000000 |
| GRN | 0.0000000 | 1.3784246 | 0.2818176 | 0.0000000 | 1.1878101 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 |
| TNF | 0.0000000 | 0.0000000 | 0.0318015 | 0.1191403 | 0.0106558 | 0.1381279 | 0.0000000 | 0.0000000 | 0.0000000 |
| LTA | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.1828448 | 0.0404983 | 0.0000000 | 0.0000000 |
| CXCL12 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.0000000 | 0.1028579 | 0.0000000 | 0.0000000 | 0.0000000 |
| CD58 | 0.0000000 | 0.0000000 | 0.1263349 | 0.1365337 | 0.0000000 | 0.1521363 | 0.0000000 | 0.1604005 | 0.0000000 |

### Build information

To keep track of the options set when running
[`build_domino()`](https://FertigLab.github.io/dominoSignal/reference/build_domino.md),
they are stored within the domino object itself. To view these options,
use the
[`dom_info()`](https://FertigLab.github.io/dominoSignal/reference/dom_info.md)
function.

``` r

dom_info(dom)
#> $create
#> [1] TRUE
#> 
#> $build
#> [1] TRUE
#> 
#> $build_variables
#>     max_tf_per_clust          min_tf_pval       max_rec_per_tf 
#>               25.000                0.001               25.000 
#> rec_tf_cor_threshold   min_rec_percentage 
#>                0.250                0.100
```

## Continued Development

Since dominoSignal is a package still being developed, there are new
functions and features that will be implemented in future versions. In
the meantime, we have put together further information on
[plotting](https://fertiglab.github.io/dominoSignal/articles/plotting_vignette)
and an example analysis can be viewed on our [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page. Additionally, if you find any bugs, have further questions, or
want to share an idea, please let us know
[here](https://github.com/FertigLab/dominoSignal/issues).

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
#> [1] dominoSignal_1.0.5
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.2.3             httr2_1.2.1           formatR_1.14         
#>   [4] biomaRt_2.66.0        rlang_1.1.6           magrittr_2.0.4       
#>   [7] clue_0.3-66           GetoptLong_1.1.0      matrixStats_1.5.0    
#>  [10] compiler_4.5.2        RSQLite_2.4.5         png_0.1-8            
#>  [13] systemfonts_1.3.1     vctrs_0.6.5           stringr_1.6.0        
#>  [16] pkgconfig_2.0.3       shape_1.4.6.1         crayon_1.5.3         
#>  [19] fastmap_1.2.0         backports_1.5.0       dbplyr_2.5.1         
#>  [22] XVector_0.50.0        rmarkdown_2.30        ragg_1.5.0           
#>  [25] purrr_1.2.0           bit_4.6.0             xfun_0.54            
#>  [28] cachem_1.1.0          jsonlite_2.0.0        progress_1.2.3       
#>  [31] blob_1.2.4            broom_1.0.11          parallel_4.5.2       
#>  [34] prettyunits_1.2.0     cluster_2.1.8.1       R6_2.6.1             
#>  [37] bslib_0.9.0           stringi_1.8.7         RColorBrewer_1.1-3   
#>  [40] car_3.1-3             jquerylib_0.1.4       Rcpp_1.1.0           
#>  [43] Seqinfo_1.0.0         iterators_1.0.14      knitr_1.50           
#>  [46] IRanges_2.44.0        Matrix_1.7-4          igraph_2.2.1         
#>  [49] tidyselect_1.2.1      dichromat_2.0-0.1     abind_1.4-8          
#>  [52] yaml_2.3.11           doParallel_1.0.17     codetools_0.2-20     
#>  [55] curl_7.0.0            lattice_0.22-7        tibble_3.3.0         
#>  [58] plyr_1.8.9            Biobase_2.70.0        withr_3.0.2          
#>  [61] KEGGREST_1.50.0       S7_0.2.1              evaluate_1.0.5       
#>  [64] desc_1.4.3            BiocFileCache_3.0.0   circlize_0.4.16      
#>  [67] Biostrings_2.78.0     pillar_1.11.1         ggpubr_0.6.2         
#>  [70] filelock_1.0.3        carData_3.0-5         foreach_1.5.2        
#>  [73] stats4_4.5.2          generics_0.1.4        S4Vectors_0.48.0     
#>  [76] hms_1.1.4             ggplot2_4.0.1         scales_1.4.0         
#>  [79] glue_1.8.0            tools_4.5.2           ggsignif_0.6.4       
#>  [82] fs_1.6.6              grid_4.5.2            tidyr_1.3.1          
#>  [85] AnnotationDbi_1.72.0  colorspace_2.1-2      Formula_1.2-5        
#>  [88] cli_3.6.5             rappdirs_0.3.3        textshaping_1.0.4    
#>  [91] ComplexHeatmap_2.26.0 dplyr_1.1.4           gtable_0.3.6         
#>  [94] rstatix_0.7.3         sass_0.4.10           digest_0.6.39        
#>  [97] BiocGenerics_0.56.0   rjson_0.2.23          htmlwidgets_1.6.4    
#> [100] farver_2.1.2          memoise_2.0.1         htmltools_0.5.9      
#> [103] pkgdown_2.2.0         lifecycle_1.0.4       httr_1.4.7           
#> [106] GlobalOptions_0.1.3   bit64_4.6.0-1
```
