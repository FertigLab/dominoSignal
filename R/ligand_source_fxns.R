
#' Find Outgoing Ligands for Each Cluster in a Domino Object
#' 
#' Generates matrices of ligand expression by each cluster and lists of each ligand expressed by a cluster over a specified threshold.
#' 
#' @param dom A domino object
#' @param signal_threshold minimum mean scaled expression of a ligand for a cell type to be considered a sender of the ligand
#' @return A list consisting of:
#' \itemize{
#'  \item{'complete_signaling'} : A matrix of mean ligand expression by each cluster
#'  \item{'binary_signaling'} : A matrix ligands by clusters where a value of 1 represents that the ligand is expressed above the signal_threshold value
#'  \iten{'outgoing_ligands'} : A list of clusters containing each ligand expressed by the cluster above the signal_threshold value
#' }
#' 
#' @examples
#' example(build_domino)
#' find_outgoing_ligands(dom = pbmc_dom_built_tiny)
#' 

find_outgoing_ligands <- function(dom, signal_threshold = 0) {
  signaling_ls <- dom@cl_signaling_matrices
  ligands <- unique(unlist(lapply(
    signaling_ls, rownames
  )))
  complete_signaling <- do.call(rbind, signaling_ls)[ligands,]
  colnames(complete_signaling) <- gsub("^L_", "", colnames(complete_signaling))
  comp_signaling <- complete_signaling
  complete_signaling[,] <- ifelse(complete_signaling > signal_threshold, 1, 0)
  binary_signaling <- complete_signaling
  dimnames(binary_signaling) <- dimnames(complete_signaling)
  senders <- colnames(binary_signaling)
  names(senders) <- senders
  outgoing_ligands <- lapply(
    senders, FUN = function(cl) {
      signal_col <- binary_signaling[, cl]
      all_lig <- names(signal_col)
      signal_logic <- signal_col > 0
      send_lig <- all_lig[signal_logic]
      return(send_lig)
    }
  )
  res <- list(
    "complete_signaling" = comp_signaling,
    "binary_signaling" = binary_signaling,
    "outgoing_ligands" = outgoing_ligands
  )
  return(res)
}