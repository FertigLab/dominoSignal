# Introducing dominoSignal: Improved Inference of Cell Signaling from Single Cell RNA Sequencing Data <a href="https://fertiglab.github.io/dominoSignal/"><img src="man/figures/logo.svg" align="right" height="138" alt="dominoSignal logo" /></a>

dominoSignal is an updated version of the original [domino](https://github.com/Elisseeff-Lab/domino) R package published in Nature Biomedical Engineering in [Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics](https://doi.org/10.1038/s41551-021-00770-5). dominoSignal is a tool for analysis of intra- and intercellular signaling in single cell RNA sequencing data based on transcription factor activation and receptor and ligand linkages between clusters.

## Installation

dominoSignal is undergoing active development to improve analysis capabilities and interpretability, so the codebase is subject to change as new features and fixes are implemented. The current version of dominoSignal can be found on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/dominoSignal.html) though the package is still undergoing development (see our [changelog](https://fertiglab.github.io/dominoSignal/news/index.html) for more information on changes). This version can be installed through Bioconductor.

```r
if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("dominoSignal")
```

The development version can be installed from the devel branch of Bioconductor.

```r
if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}

# The following initializes usage of Bioc devel
BiocManager::install(version = "devel")
BiocManager::install("dominoSignal")
```

## Usage Overview

Here is an overview of how dominoSignal might be used in analysis of a single cell RNA sequencing data set:

1. Transcription factor activation scores are calculated (we recommend using [pySCENIC](https://pyscenic.readthedocs.io/en/latest/), but other methods can be used as well)
2. A ligand-receptor database is used to map linkages between ligands and receptors (we recommend using [CellPhoneDB](https://www.cellphonedb.org/), but other methods can be used as well).
3. A domino object is created using counts, z-scored counts, clustering information, and the data from steps 1 and 2.
4. Parameters such as the maximum number of transcription factors and receptors or the minimum correlation threshold (among others) are used to make a cell communication network
5. Communication networks can be extracted from within the domino object or visualized using a variety of plotting functions

Please see [our website](https://fertiglab.github.io/dominoSignal/) for tutorials on all of these steps, from downloading and running [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) in the [SCENIC tutorial](https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette.html) to building and visualizing domino results on the [Getting Started page](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html). Other articles include [further details on plotting functions](https://fertiglab.github.io/dominoSignal/articles/plotting_vignette.html) and [the structure of the domino object](https://fertiglab.github.io/dominoSignal/articles/domino_object_vignette.html).

## Citation

If you use the development version of our package in your analysis, please cite us:

> Cherry C, Maestas DR, Han J, Andorko JI, Cahan P, Fertig EJ, Garmire LX, Elisseeff JH. Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics. Nat Biomed Eng. 2021 Oct;5(10):1228-1238. doi: 10.1038/s41551-021-00770-5. Epub 2021 Aug 2. PMID: 34341534; PMCID: PMC9894531.
> Cherry C, Mitchell J, Nagaraj S, Krishnan K, Lvovs D, Fertig E, Elisseeff J (2026). dominoSignal: Cell Communication Analysis for Single Cell RNA Sequencing. R package version 1.7.0.

## Contact Us

If you find any bugs or have questions, please [open an issue here](https://github.com/FertigLab/dominoSignal/issues).
