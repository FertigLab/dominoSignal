
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
#'  \item{'outgoing_ligands'} : A list of clusters containing each ligand expressed by the cluster above the signal_threshold value
#' }
#' @export
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

#' Query Ligand Senders
#' 
#' Identify clusters that are senders of ligands received by a specified cluster
#' 
#' @param dom A domino object
#' @param outgoing_list list output from [find_outgoing_ligands()]
#' \itemize{
#'  \item{'complete_signaling'} : A matrix of mean ligand expression by each cluster
#'  \item{'binary_signaling'} : A matrix ligands by clusters where a value of 1 represents that the ligand is expressed above the signal_threshold value
#'  \item{'outgoing_ligands'} : A list of clusters containing each ligand expressed by the cluster above the signal_threshold value
#' }
#' @param receiver_cell name of the receiving cluster whose receptors will be used 
#' @return list of incoming ligands to the recipient cell type organized by their receptor
#' @export
#' 
#' @examples
#' example(build_domino)
#' outgoing_res <- find_outgoing_ligands(dom = pbmc_dom_built_tiny)
#' query_ligand_senders(
#'   dom = pbmc_dom_built_tiny,
#'   outgoing_list = outgoing_res,
#'   receiver_cell = "CD14_monocyte"
#' )
#' 

query_ligand_senders <- function(dom, outgoing_list, receiver_cell) {
  clusters <- names(dom@linkages$clust_rec)
  # find active receptors
  cl_rec <- dom@linkages$clust_rec[[receiver_cell]]
  rec_lig_ls <- dom@linkages$rec_lig[cl_rec]
  sender_ls <- lapply(
    rec_lig_ls, 
    FUN = function(ligs) {
      cl_sender_ls <- lapply(
        clusters, FUN = function(cl) {
          return(as.numeric(ligs %in% outgoing_list$outgoing_ligands[[cl]]))
        }
      )
      cl_sender <- do.call(cbind, cl_sender_ls)
      dimnames(cl_sender) <- list(
        ligs,
        clusters
      )
      return(cl_sender)
    }
  )
  return(sender_ls)
}

#' Create a vector of intercellular linkages
#' 
#' Rephrase an incoming ligand query into a vector of intercellular linkages in the format of:
#' '\["receiver"\]:"receptor" <- "ligand":\["sender"\]'
#' 
#' @param ligand_query list of incoming ligands to the recipient cell type organized by their receptor
#' @param receiver_cell name of the receiving cluster
#' @param receptor name of the receptor receiving signals
#' @return vector of intercellular linkages for each sender cell and active ligand
#' @export
#' 
#' @examples
#' example(build_domino)
#' LQ <- query_ligand_senders(
#'   dom = pbmc_dom_built_tiny,
#'   outgoing_list = find_outgoing_ligands(dom = pbmc_dom_built_tiny),
#'   receiver_cell = "CD14_monocyte"
#' )
#' senders_as_vector(
#'   ligand_query = LQ,
#'   receiver_cell = "CD14_monocyte",
#'   receptor = "CXCR3"
#' )
#' 
#' 
senders_as_vector <- function(ligand_query, receiver_cell, receptor) {
  ligand_senders <- ligand_query[[receptor]]
  ligands <- rownames(ligand_senders)
  senders <- colnames(ligand_senders)
  linkage <- c()
  for(i in seq_along(ligands)) {
    lig <- ligands[i]
    for(j in seq_along(senders)) {
      cl <- senders[j]
      if(ligand_senders[i,j] == 1) {
        linkage <- c(
          linkage,
          paste0("[", receiver_cell, "]:", receptor, " <- ", lig, ":[", cl, "]")
        )
      }
    }
  }
  return(linkage)
}

#' Add intercellular linkages to a domino object
#' 
#' @param dom A domino object
#' @param signal_threshold minimum mean scaled expression of a ligand for a cell type to be considered a sender of the ligand
#' 
#' @return a domino object including intercellular linkages within the linkages slot as "rec_lig_cl"
#' @export
#' 
#' @examples
#' example(build_domino)
#' dom <- add_intercellular_linkages(dom = pbmc_dom_built_tiny, signal_threshold = 0)
#' slot(dom, "linkages")$rec_lig_cl
#' 

add_intercellular_linkages <- function(dom, signal_threshold = 0){
  rec_lig_cl <- list()
  cell_types <- levels(dom@clusters)
  for(j in seq_along(cell_types)) {
    cl <- cell_types[j]
    sender_query <- query_ligand_senders(
      dom, 
      outgoing_list = find_outgoing_ligands(dom), 
      receiver_cell = cl
    )
    cl_recs <- names(sender_query)
    names(cl_recs) <- cl_recs
    sender_vec <- unlist(
      lapply(cl_recs, FUN = function(x) {senders_as_vector(sender_query, receiver_cell = cl, receptor = x)}),
      use.names = FALSE
    )
    rec_lig_cl[[cl]] <- sender_vec
  }
  dom@linkages$rec_lig_cl <- rec_lig_cl
  return(dom)
}


