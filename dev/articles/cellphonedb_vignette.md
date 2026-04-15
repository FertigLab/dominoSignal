# Using the CellPhoneDB Database

The dominoSignal analysis pipeline requires the use of a reference
receptor - ligand database to identify potential receptor - ligand
interactions. We recommend the use of CellPhoneDB. Here is a very brief
tutorial on how to download the requisite files from CellPhoneDB v4.0.0.

## File Downloads for CellPhoneDB:

Database files from [CellPhoneDB
v4.0.0](https://github.com/ventolab/cellphonedb-data/releases/tag/v4.0.0)
for human scRNAseq data can be installed from a public Github repository
from the Tiechmann Group that developed CellPhoneDB.

``` r
# URL for desired version of CellPhoneDB
cellphone_url <- "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz"

# download compressed database
cellphone_tar <- paste0(temp_dir, "/cellphoneDB_v4.tar.gz")
download.file(url = cellphone_url, destfile = cellphone_tar)

# move contents of the compressed file to a new directory
untar(tarfile = cellphone_tar, exdir = cellphone_dir)
cellphone_data <- paste0(cellphone_dir, "/cellphonedb-data-4.0.0/data")
```

To facilitate the use of these files in the format used in dominoSignal,
we include a helper function,
[`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_rl_map_cellphonedb.md),
that automatically parses files from the CellPhoneDB database to arrive
at the rl_map format. For more information on how to use these files in
the dominoSignal pipeline, please see our [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html)
page. To learn how to use SCENIC for TF activation scoring, please see
our [SCENIC for TF Activation
Scoring](https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette)
tutorial.

Vignette Build Information

Date last built and session information:

``` r
Sys.Date()
#> [1] "2026-04-15"
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       formatR_1.14      otel_0.2.0       
#> [25] jsonlite_2.0.0    rlang_1.2.0       fs_2.0.1          htmlwidgets_1.6.4
```
