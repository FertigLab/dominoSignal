---
output: github_document
---

<!-- index.md is generated from index.Rmd. Please edit that file -->



## Introducing dominoSignal <a href="https://fertiglab.github.io/dominoSignal/"><img src="man/figures/logo.svg" align="right" height="138" style="float:right; height:138px;" alt="dominoSignal logo" /></a>

dominoSignal is an updated version of the original [domino](https://github.com/Elisseeff-Lab/domino) R package published in Nature Biomedical Engineering in [Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics](https://doi.org/10.1038/s41551-021-00770-5). dominoSignal is a tool for analysis of intra- and intercellular signaling in single cell RNA sequencing data based on transcription factor activation and receptor and ligand linkages.

## Installation

This version is currently hosted on the [FertigLab GitHub](https://github.com/FertigLab) on the [dominoSignal repository](https://github.com/FertigLab/dominoSignal) forked from the primary repository hosted on the [Elisseeff-Lab GitHub](https://github.com/Elisseeff-Lab/domino), and can be installed using the remotes package. This is the current stable version during these active updates for reproducible usage (see our [changelog](news/index.html) for more information on changes).


```r
if (!require(remotes)) {
    install.packages("remotes")
}
remotes::install_github("FertigLab/dominoSignal")
```

## Usage Overview

Here is an overview of how dominoSignal might be used in analysis of a single cell RNA sequencing data set:

1. Transcription factor activation scores are calculated (we recommend using [pySCENIC](https://pyscenic.readthedocs.io/en/latest/), but other methods can be used as well). For more information on how to use SCENIC, please see our [Using SCENIC for TF Activation](https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette) page.
2. A ligand-receptor database is used to map linkages between ligands and receptors (we recommend using [cellphoneDB](https://www.cellphonedb.org/), but other methods can be used as well). For information on downloading the necessary files for cellphoneDB, please see our [Using the cellphoneDB Database](https://fertiglab.github.io/dominoSignal/articles/cellphonedb_vignette) page.
3. A domino object is created using counts, z-scored counts, clustering information, and the data from steps 1 and 2.
4. Parameters such as the maximum number of transcription factors and receptors or the minimum correlation threshold (among others) are used to make a cell communication network
5. Communication networks can be extracted from within the domino object or visualized using a variety of plotting functions

Please see the [Getting Started](https://fertiglab.github.io/dominoSignal) page for an example analysis that includes all of these steps in a dominoSignal analysis in detail, from creating a domino object to parameters for building the network to visualizing domino results. Other articles include further details on [plotting functions](https://fertiglab.github.io/dominoSignal/articles/plotting_vignette) and the structure of the [domino object](https://fertiglab.github.io/dominoSignal/articles/domino_object_vignette).
