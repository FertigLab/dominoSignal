#' @import plyr
#' @import methods
#'
NULL

#' Renames clusters in a domino object
#'
#' This function renames the clusters used to build a domino object
#'
#' @param dom a domino object to rename clusters in
#' @param clust_conv named vector of conversions from old to new clusters. Values are taken as new clusters IDs and names as old cluster IDs.
#' @param warning logical. If TRUE, will warn if a cluster is not found in the conversion table. Default is FALSE.
#' @return A domino object with clusters renamed in all applicable slots.
#' @export
#' @examples 
#' example(build_domino, echo = FALSE)
#' new_clust <- c("CD8_T_cell" = "CD8+ T Cells",
#'  "CD14_monocyte" = "CD14+ Monocytes", "B_cell" = "B Cells")
#' pbmc_dom_built_tiny <- rename_clusters(pbmc_dom_built_tiny, new_clust)
#'
rename_clusters <- function(dom, clust_conv, warning = FALSE) {
    if (is.null(dom@clusters)) {
        stop("There are no clusters in this domino object")
    }
    if (dom@misc$create) {
        dom@clusters <- plyr::revalue(dom@clusters, clust_conv, warn_missing = warning)
        colnames(dom@clust_de) <- plyr::revalue(colnames(dom@clust_de), clust_conv, warn_missing = warning)
        names(colnames(dom@clust_de)) <- c()
        colnames(dom@misc$cl_rec_percent) <- plyr::revalue(colnames(dom@misc$cl_rec_percent),
            clust_conv,
            warn_missing = warning
        )
    }
    if (dom@misc$build) {
        names(dom@linkages$clust_tf) <- plyr::revalue(names(dom@linkages$clust_tf), clust_conv, warn_missing = warning)
        names(dom@linkages$clust_rec) <- plyr::revalue(names(dom@linkages$clust_rec), clust_conv, warn_missing = warning)
        names(dom@linkages$clust_incoming_lig) <- plyr::revalue(names(dom@linkages$clust_incoming_lig),
            clust_conv,
            warn_missing = warning
        )
        names(dom@linkages$clust_tf_rec) <- plyr::revalue(names(dom@linkages$clust_tf_rec),
            clust_conv,
            warn_missing = warning
        )
        sig_ligands <- colnames(dom@signaling)
        sig_rec <- rownames(dom@signaling)
        # Remove L_ prefix
        sig_ligand_clust <- gsub("^L_", "", sig_ligands)
        # Remove R_ prefix
        sig_rec_clust <- gsub("^R_", "", sig_rec)
        new_lig_clust <- plyr::revalue(sig_ligand_clust, clust_conv, warn_missing = warning)
        new_rec_clust <- plyr::revalue(sig_rec_clust, clust_conv, warn_missing = warning)
        colnames(dom@signaling) <- paste0("L_", new_lig_clust)
        rownames(dom@signaling) <- paste0("R_", new_rec_clust)
        names(dom@cl_signaling_matrices) <- plyr::revalue(names(dom@cl_signaling_matrices),
            clust_conv,
            warn_missing = warning
        )
        for (cl in names(dom@cl_signaling_matrices)) {
            cl_sig_ligands <- colnames(dom@cl_signaling_matrices[[cl]])
            # Remove L_ prefix
            cl_sig_lig_clust <- gsub("^L_", "", cl_sig_ligands)
            cl_sig_lig_new <- plyr::revalue(cl_sig_lig_clust, clust_conv, warn_missing = warning)
            colnames(dom@cl_signaling_matrices[[cl]]) <- paste0("L_", cl_sig_lig_new)
        }
    }
    return(dom)
}

#' Adds a column to the RL signaling data frame.
#'
#' This function adds a column to the internal rl 'map' used to map all
#' receptor and receptor complexes to all ligand and ligand complexes.
#'
#' @param map RL signaling data frame.
#' @param map_ref Name of column to match new data to
#' @param conv Data frame matching current data in map to new data.
#' @param new_name Name of new column to be created in RL map
#' @return An updated RL signaling data frame
#' @export
#' @examples 
#' example(create_rl_map_cellphonedb, echo = FALSE)
#' lr_name <- data.frame("abbrev" = c("L", "R"), "full" = c("Ligand", "Receptor"))
#' rl_map_expanded <- add_rl_column(map = rl_map_tiny, map_ref = "type_A",
#' conv = lr_name, new_name = "type_A_full")
#' 
add_rl_column <- function(map, map_ref, conv, new_name) {
  map_in_ref <- match(map[[map_ref]], conv[, 1])
  not_in_ref <- which(is.na(map_in_ref))
  if (length(not_in_ref > 0)) {
    not_in_ref_map <- cbind.data.frame(map[not_in_ref, ], as.character(NA), stringsAsFactors = FALSE)
    colnames(not_in_ref_map)[ncol(not_in_ref_map)] <- new_name
    rownames(not_in_ref_map) <- c()
  } else {
    not_in_ref_map <- c()
  }
  new_map <- c()
  for (r_id in seq_len(nrow(map))) {
    row <- map[r_id, ]
    conv_ids <- which(conv[, 1] == row[[map_ref]])
    for (id in conv_ids) {
      new_row <- c(as.matrix(row), conv[id, 2])
      new_map <- rbind(new_map, new_row)
    }
  }
  rownames(new_map) <- c()
  colnames(new_map) <- c(colnames(map), new_name)
  new_map <- rbind.data.frame(new_map, not_in_ref_map, stringsAsFactors = FALSE)
  new_map <- data.frame(new_map, stringsAsFactors = FALSE)
}