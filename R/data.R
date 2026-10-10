#' SCENIC example subset
#'
#' A subset of pySCENIC output for the PBMC example data.
#'
#' @format A list of:
#' \describe{
#'  \item{auc_tiny}{Matrix of AUC scores from pySCENIC for a subset of TFs and cells.}
#'  \item{regulons_tiny}{Data frame subset of the pySCENIC regulon table.}
#'  \item{regulon_list_tiny}{The pySCENIC regulons formatted as a named list by [create_regulon_list_scenic()].}
#' }
#'
#' @source <https://doi.org/10.5281/zenodo.10124865> (files `auc_pbmc_3k.csv`, `regulons_pbmc_3k.csv`); formatted with `data-raw/tiny_data.R`
#' @usage data("SCENIC")
"SCENIC"


#' PBMC scRNAseq data subset
#'
#' A subset of the preprocessed 10x Genomics PBMC 3k scRNAseq data.
#'
#' @format A list of:
#' \describe{
#'  \item{count_tiny}{Subset of the sparse counts matrix (`dgCMatrix`), genes x cells.}
#'  \item{zscore_tiny}{Subset of z-scored expression matrix, genes x cells.}
#'  \item{clusters_tiny}{Factor of cell-type labels named by cell barcode.}
#' }
#'
#' @source <https://doi.org/10.5281/zenodo.10124865> (file `pbmc3k_sce.rds`); formatted with `data-raw/tiny_data.R`.
#' @usage data("PBMC")
"PBMC"


#' CellPhoneDB subset
#'
#' A list of four subsets of CellPhoneDB v4.0.0 data and resulting domino rl_map.
#'
#' @format A list of:
#' \describe{
#'  \item{genes_tiny}{A subset of the CellPhoneDB `gene_input.csv` file.}
#'  \item{proteins_tiny}{A subset of the CellPhoneDB `protein_input.csv` file.}
#'  \item{complexes_tiny}{A subset of the CellPhoneDB `complex_input.csv` file.}
#'  \item{interactions_tiny}{A subset of the CellPhoneDB `interaction_input.csv` file.}
#'  \item{rl_map_tiny}{The subset formatted as a receptor-ligand map by [create_rl_map_cellphonedb()].}
#' }
#' 
#' @source <https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz>; subset with `data-raw/tiny_data.R`.
#' @usage data("CellPhoneDB")
"CellPhoneDB"

#' Example domino objects
#' 
#' A list of two domino objects built from the example data subsets ([PBMC], [SCENIC], and [CellPhoneDB]).
#' \describe{
#'  \item{dom_tiny}{A [domino()] object created using [create_domino()].}
#'  \item{built_dom_tiny}{The same object after running [build_domino()].}
#' }
#' @source Generated with `data-raw/tiny_data.R`.
#' @usage data("DominoObjects")
"DominoObjects"

#' Example linkage summary
#' 
#' A list containing a mock linkages summary across six subjects in two groups,
#'   and the result of running [test_differential_linkages()] on it.
#' @format A list of:
#' \describe{
#' \item{linkage_sum_tiny}{A [linkage_summary()] object created with mock data.}
#' \item{linkage_diff_tiny}{Results from running [test_differential_linkages()] for receptors (`linkage = "rec"`) in cluster `"C1"`, grouped by `"group"`.}
#' }
#' @source Generated with `data-raw/tiny_data.R`.
#' @usage data("LinkageSummary")
"LinkageSummary"