# Using the CellPhoneDB Database

The dominoSignal analysis pipeline requires the use of a reference
receptor-ligand database to identify potential receptor-ligand
interactions. We recommend CellPhoneDB. Here is a very brief tutorial on
how to download the requisite files from CellPhoneDB v4.0.0.

## File downloads for CellPhoneDB:

Database files from [CellPhoneDB
v4.0.0](https://github.com/ventolab/cellphonedb-data/releases/tag/v4.0.0)
for human scRNAseq data can be installed from a public GitHub repository
from the Teichmann Group, which developed CellPhoneDB.

``` r

# directories for the download
temp_dir <- tempdir()
cellphone_dir <- file.path(temp_dir, "cellphone")

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
at the rl_map format.

## Related articles

- [Get started with
  dominoSignal](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html):
  the basic workflow
- [Interacting with domino
  objects](https://fertiglab.github.io/dominoSignal/articles/domino_object_vignette.html):
  object structure and accessor functions
- [Plotting functions and
  options](https://fertiglab.github.io/dominoSignal/articles/plotting_vignette.html):
  visualization options
- [Differential signaling
  workflow](https://fertiglab.github.io/dominoSignal/articles/differential_signaling.html):
  comparing signaling across subjects
- [SCENIC for TF activity
  scoring](https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette.html):
  generating TF activity inputs

------------------------------------------------------------------------

To report a bug, ask a question, or suggest a feature, please [open an
issue](https://github.com/FertigLab/dominoSignal/issues).

Vignette Build Information

Date last built and session information:

``` r

Sys.Date()
#> [1] "2026-11-24"
sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
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
#>  [5] xfun_0.62         cachem_1.1.0      knitr_1.52        htmltools_0.5.9  
#>  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         sass_0.5.0       
#> [13] pkgdown_2.2.1     textshaping_1.1.0 jquerylib_0.1.4   systemfonts_1.3.4
#> [17] compiler_4.6.2    tools_4.6.2       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       formatR_1.14      otel_0.3.0       
#> [25] jsonlite_2.0.0    rlang_1.3.0       fs_2.1.0
```
