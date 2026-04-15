# SCENIC for TF Activation Scoring

[SCENIC](https://scenic.aertslab.org/) is our preferred software for TF
activity scoring. We recommend using the Python implementation
([pySCENIC](https://pyscenic.readthedocs.io/en/latest/)) as it is faster
than the original R implementation. Be aware, the use of
[SCENIC](https://scenic.aertslab.org/) is often the slowest and most
memory intensive step of this analysis pipeline.
[SCENIC](https://scenic.aertslab.org/) should be run on computing
resources with access to multi-core processing and large amounts of
memory.

## File Downloads for SCENIC

For this tutorial,
[pySCENIC](https://pyscenic.readthedocs.io/en/latest/) is used as the
method for TF activity inference. This requires downloading some files
prior to use. Here, we show how to download these files and the way we
run pySCENIC to generate the necessary files for use in dominoSignal
analyis.

A singularity image of [SCENIC](https://scenic.aertslab.org/) v0.12.1
can be installed from the DockerHub image as a source.
[SCENIC](https://scenic.aertslab.org/) requires a list of TFs, motif
annotations, and cisTarget motifs which are all available from the
authors of [SCENIC](https://scenic.aertslab.org/) for human (HGNC),
mouse (MGI), and fly. The following will download everything necessary
for an analysis of a data set with HGNC gene labels for the hg38 genome.

### Reference files

First we set up a temporary directory to store the files.

``` r
temp_dir <- tempdir()
```

Then we build the singularity image and download reference files to the
temporary directory.

``` bash
SCENIC_DIR="'temp_dir'/scenic"
mkdir ${SCENIC_DIR}

# Build singularity image
singularity build "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" docker://aertslab/pyscenic:0.12.1

# Matrix containing motifs as rows and genes as columns and ranking position for each gene and motif (based on CRM scores) as values. URLs provided link to v2 feather files required for the 0.12.1 version of pySENIC.

curl "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  -o "${SCENIC_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

curl "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  -o "${SCENIC_DIR}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

# List of genes encoding TFs in the HG38 reference genome
curl "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt" \
  -o "${SCENIC_DIR}/allTFs_hg38.txt"

# Motif annotations based on the 2017 cisTarget motif collection. Use these files if you are using the mc9nr databases.
curl "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
  -o "${SCENIC_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
```

### Preparing Preprocessed Data:

dominoSignal was designed to be compatible with single cell workflows,
accepting gene by cell matrices for counts and scaled counts. However,
[pySCENIC](https://pyscenic.readthedocs.io/en/latest/) is used for TF
activation inference, and it requires as input a cell by gene matrix.
This is the opposite orientation of many R based single cell packages
(such as [Seurat](https://satijalab.org/seurat) or
[SingleCellExperiment](https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment))
but is the default for Python based tools such as
[scanpy](https://scanpy.readthedocs.io/en/stable/). The RNA counts
matrix can be saved as a tab-seperated value (.tsv) or comma-sperated
value (.csv) file after transposing the matrix to cell by gene
orientation. Below is an example of extracting the counts from a
[SingleCellExperiment](https://doi.org/doi:10.18129/B9.bioc.SingleCellExperiment)
object and saving it as a .tsv file.

``` r
pbmc_counts <- assay(pbmc, "counts")
write.table(t(as.matrix(pbmc_counts)), paste0(input_dir, "/pbmc3k_counts.tsv"), sep = "\t",
    col.names = NA)
```

tsv and csv files are very inefficient for storing large matrices. As an
alternative, the counts matrix can be saved as a loom file and directly
passed to [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) in this
format. Generating a loom file in R requires use of the
[loomR](https://github.com/mojaveazure/loomR) package. The package
automatically handles the transposition of the counts matrix to cell by
gene orientation as well. We recommend this approach to passing a counts
matrix to pySCENIC. However, users should be aware that
[loomR](https://github.com/mojaveazure/loomR) is not hosted on CRAN or
Bioconductor at the time of this vignette’s creation.

``` r
library(loomR)

# save loom counts matrix
pbmc_counts <- assay(pbmc, "counts")
pbmc_loom <- loomR::create(filename = paste0(input_dir, "/pbmc3k_counts.loom"), data = pbmc_counts)
# Depending on the version of loomR, you may need to manually add CellId and
# Gene attributes:
pbmc_loom$add.row.attribute(list(Gene = rownames(pbmc)))
pbmc_loom$add.col.attribute(list(CellID = colnames(pbmc)))

# connection to the loom file must be closed to complete file creation
pbmc_loom$close_all()
```

We will use the `pbmc3k\_counts.loom` file in this tutorial.

## Running SCENIC

[pySCENIC](https://pyscenic.readthedocs.io/en/latest/) is initiated
using bash scripting in the terminal. The analysis consists of three
steps to score genes for TF motif enrichment, construct TF regulons
consisting of genes targeted by the TFs, and arrive at AUC scores for
enrichment of regulon gene transcription within each cell.

### grn: construct TF-modules

Co-expression modules are used to quantify gene-TF adjacencies.

``` bash
singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic grn \
    "pbmc3k_counts.loom" \
    "${SCENIC_DIR}/allTFs_hg38.txt" \
    -o "${SCENIC_DIR}/pbmc_adj_3k.tsv" \
    --num_workers 6 \
    --seed 123

# Arguments:
    # path to the loom file
    # list of TFs
    # output directory for the adjacency matrix
    # number of CPUs to use if multi-core processing is available
    # specify a random seed for reproducibility
```

### ctx: construct TF regulons with pruning based on TF motif enrichment

The rankings of genes based on enrichment of TF motifs to the
transcription start site (TSS) are considered in the construction of
regulons, where target genes in the TF-modules are removed if they lack
motifs proximal to the TSS where the TF may bind.

``` bash
singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic ctx \
    "${SCENIC_DIR}/pbmc_adj.tsv" \
    "${SCENIC_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "${SCENIC_DIR}/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "${SCENIC_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname "pbmc3k_counts.loom" \
    --mode "dask_multiprocessing" \
    --output "${SCENIC_DIR}/regulons_pbmc_3k.csv" \
    --num_workers 1                                                          

# Arguments:
    # adjacency matrix output from grn
    # target rankings of motif enrichment within 10 kb of TSS
    # target rankings of motif enrichment within 500 bp upstream and 100 bp downstream of TSS
    # TF motif features
    # counts matrix loom file
    # enable multi-core processing
    # output file of learned TF regulons
    # number of CPU cores
```

### aucell: calculate TF activity scores

Enrichment of a regulon is measured as the **A**rea **U**nder the
recovery **C**urve (AUC) of the genes that define this regulon.

``` bash
singularity exec "${SCENIC_DIR}/aertslab-pyscenic-0.12.1.sif" pyscenic aucell \
    "pbmc3k_counts.loom" \
    "${SCENIC_DIR}/regulons_pbmc_3k.csv" \
    -o "${SCENIC_DIR}/auc_pbmc_3k.csv"
    
# Arguments:
    # counts matrix loom file
    # regulon table output from ctx
    # cell x TF matrix of TF enrichment AUC values
```

The csv files created in this analysis can be used as the TF activation
inputs in a dominoSignal analysis. Please see our [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html)
page for more information on how to perform an analysis with
dominoSignal. For a brief tutorial on the download of the necessary
files for using CellPhoneDB, please see the [Using the CellPhoneDB
database](https://fertiglab.github.io/dominoSignal/articles/cellphonedb_vignette.html)
vignette.

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
