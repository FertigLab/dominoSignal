#' SCENIC AUC subset
#'
#' A subset of SCENIC AUCs as applied to PBMC data.
#'
#' @format A list of:
#' \describe{
#'  \item{auc_tiny}{A subset of SCENIC AUCs}
#'  \item{regulons_tiny}{A subset of SCENIC regulons}
#'  \item{regulon_list_tiny}{A subset of SCENIC regulons formatted as a regulon_list}
#' }
#'
#' @source <https://zenodo.org/records/10951634/files>
#' @usage data("SCENIC")
"SCENIC"


#' PBMC RNAseq data subset
#'
#' A subset of the results of PBMC RNA-seq data.
#'
#' @format A list of::
#' \describe{
#'  \item{count_tiny}{A subset of PBMC RNA-seq data: counts assay}
#'  \item{zscore_tiny}{A subset of PBMC RNA-seq data: zscore assay}
#'  \item{clusters_tiny}{A subset of PBMC RNA-seq data: clusters as defined by cell_type}
#' }
#'
#' @source <https://zenodo.org/records/10951634/files/pbmc3k_sce.rds>
#' @usage data("PBMC")
"PBMC"


#' CellPhoneDB subset
#'
#' A list of four subsets of CellPhoneDB data.
#'
#'
#' @format A list of:
#' \describe{
#'  \item{genes_tiny}{A subet of CellPhoneDB gene_input.csv}
#'  \item{proteins_tiny}{A subset of CellPhoneDB protein_input.csv}
#'  \item{complexes_tiny}{A subset of CellPhoneDB complex_input.csv}
#'  \item{interactions_tiny}{A subset of CellPhoneDB interaction_input.csv}
#'  \item{rl_map_tiny}{A subset of CellPhoneDB data formatted as a domino rl_map}
#' }
#' 
#' @source <https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz>
#' @usage data("CellPhoneDB")
"CellPhoneDB"

#' Example domino objects
#' 
#' A list of two domino objects, one from [create_domino()] and one from [build_domino()] outputs.
#' \describe{
#'  \item{dom_tiny}{A domino object created using [create_domino()] with the tiny datasets.}
#'  \item{dom_tiny_built}{A domino object created using [build_domino()] with the tiny datasets.}
#' }
#' @usage data("DominoObjects")
"DominoObjects"

#' Example linkage summary
#' 
#' A mock linkage summary object with made up data to show the structure of the object class
#' 
#' @usage data("LinkageSummary")
"LinkageSummary"