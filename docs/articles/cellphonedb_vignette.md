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
[`create_rl_map_cellphonedb()`](https://FertigLab.github.io/dominoSignal/reference/create_rl_map_cellphonedb.md),
that automatically parses files from the CellPhoneDB database to arrive
at the rl_map format. For more information on how to use these files in
the dominoSignal pipeline, please see our [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal)
page. To learn how to use SCENIC for TF activation scoring, please see
our [SCENIC for TF Activation
Scoring](https://fertiglab.github.io/dominoSignal/articles/articles/tf_scenic_vignette)
tutorial.

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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.54         cachem_1.1.0      knitr_1.50        htmltools_0.5.9  
#>  [9] rmarkdown_2.30    lifecycle_1.0.4   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
#> [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        bslib_0.9.0      
#> [21] evaluate_1.0.5    yaml_2.3.11       formatR_1.14      jsonlite_2.0.0   
#> [25] rlang_1.1.6       fs_1.6.6          htmlwidgets_1.6.4
```
