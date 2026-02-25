## Introducing dominoSignal [![dominoSignal logo](logo.svg)](https://fertiglab.github.io/dominoSignal/)

dominoSignal is an updated version of the original
[domino](https://github.com/Elisseeff-Lab/domino) R package published in
Nature Biomedical Engineering in [Computational reconstruction of the
signalling networks surrounding implanted biomaterials from single-cell
transcriptomics](https://doi.org/10.1038/s41551-021-00770-5).
dominoSignal is a tool for analysis of intra- and intercellular
signaling in single cell RNA sequencing data based on transcription
factor activation and receptor and ligand linkages.

## Installation

dominoSignal is undergoing active development to improve analysis
capabilities and interpretability, so the codebase is subject to change
as new features and fixes are implemented. The current version of
dominoSignal can be found on
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/dominoSignal.html)
though the package is still undergoing development (see our
[changelog](https://FertigLab.github.io/dominoSignal/news/index.md) for
more information on changes). This version can be installed through
Bioconductor.

``` r

if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("dominoSignal")
```

The development version is currently hosted on the [FertigLab
GitHub](https://github.com/FertigLab) on the [dominoSignal GitHub
repository](https://github.com/FertigLab/dominoSignal), and can be
installed using the remotes package.

``` r

if (!require(remotes)) {
    install.packages("remotes")
}
remotes::install_github("FertigLab/dominoSignal")
```

## Usage Overview

Here is an overview of how dominoSignal might be used in analysis of a
single cell RNA sequencing data set:

1.  Transcription factor activation scores are calculated (we recommend
    using [pySCENIC](https://pyscenic.readthedocs.io/en/latest/), but
    other methods can be used as well). For more information on how to
    use SCENIC, please see our [Using SCENIC for TF
    Activation](https://FertigLab.github.io/dominoSignal/vignette(%22articles/scenic_vignette%22))
    page.
2.  A ligand-receptor database is used to map linkages between ligands
    and receptors (we recommend using
    [cellphoneDB](https://www.cellphonedb.org/), but other methods can
    be used as well). For information on downloading the necessary files
    for cellphoneDB, please see our [Using the cellphoneDB
    Database](https://FertigLab.github.io/dominoSignal/vignette(%22articles/cellphonedb_vignette%22))
    page.
3.  A domino object is created using counts, z-scored counts, clustering
    information, and the data from steps 1 and 2.
4.  Parameters such as the maximum number of transcription factors and
    receptors or the minimum correlation threshold (among others) are
    used to make a cell communication network
5.  Communication networks can be extracted from within the domino
    object or visualized using a variety of plotting functions

Please see the [Getting
Started](https://FertigLab.github.io/dominoSignal/vignette(%22dominoSignal%22))
page for an example analysis that includes all of these steps in a
dominoSignal analysis in detail, from creating a domino object to
parameters for building the network to visualizing domino results. Other
articles include further details on [plotting
functions](https://FertigLab.github.io/dominoSignal/vignette(%22plotting_vignette%22))
and the structure of the [domino
object](https://FertigLab.github.io/dominoSignal/vignette(%22domino_object_vignette%22)).
