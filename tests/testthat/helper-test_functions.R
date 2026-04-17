#' Make a domino object with duplicate ligand aliases for testing
#' 
#' @param dom A domino object to modify
#' @param alias_names A vector of ligand alias names to add
#' @param target_gene The gene name that the aliases should point to
#' @param receptor_slots A vector of receptor names to which the aliases should be linked
#' @return A modified domino object with duplicate ligand aliases

dom_dup_ligand_aliases <- function(dom, alias_names, target_gene, receptor_slots) {
    test_dom <- dom
    template_row <- test_dom@misc[["rl_map"]][1, , drop = FALSE]
    alias_rows <- template_row[rep(1, length(alias_names)), , drop = FALSE]
    alias_rows$L.name <- alias_names
    alias_rows$L.gene <- rep(target_gene, length(alias_names))
    test_dom@misc[["rl_map"]] <- rbind(test_dom@misc[["rl_map"]], alias_rows)

    for (i in seq_along(alias_names)) {
        test_dom@linkages$rec_lig[[receptor_slots[i]]] <- alias_names[i]
    }

    return(test_dom)
}

#' Make a domino object with duplicate receptor aliases for testing
#' 
#' @param dom A domino object to modify
#' @param alias_names A vector of receptor alias names to add
#' @param target_gene The gene name that the aliases should point to
#' @param ligand_slots A vector of ligand names to which the aliases should be linked
#' @return A modified domino object with duplicate receptor aliases

dom_dup_receptor_aliases <- function(dom, alias_names, target_gene, ligand_slots) {
    test_dom <- dom
    template_row <- test_dom@misc[["rl_map"]][1, , drop = FALSE]
    alias_rows <- template_row[rep(1, length(alias_names)), , drop = FALSE]
    alias_rows$R.name <- alias_names
    alias_rows$R.gene <- rep(target_gene, length(alias_names))
    test_dom@misc[["rl_map"]] <- rbind(test_dom@misc[["rl_map"]], alias_rows)

    for (i in seq_along(alias_names)) {
        test_dom@linkages$rec_lig[[alias_names[i]]] <- ligand_slots[i]
    }

    return(test_dom)
}

#' Add zero expression values for a specified gene in a domino object
#' 
#' @param dom A domino object to modify
#' @param gene The gene name for which to add zero expression values
#' @return A modified domino object with zero expression values for the specified gene

dom_add_zero_exp <- function(dom, gene) {
    if (!(gene %in% rownames(dom@counts)) || !(gene %in% rownames(dom@z_scores))) {
        stop("Gene not found in domino object")
    }
    test_dom <- dom
    test_dom@counts[gene, ] <- 0
    test_dom@z_scores[gene, ] <- 0
    return(test_dom)
}

#' Add zero AUC values for a specified transcription factor in a domino object
#' 
#' @param dom A domino object to modify
#' @param tf The transcription factor name for which to add zero AUC values
#' @return A modified domino object with zero AUC values for the specified transcription factor

dom_add_zero_auc <- function(dom, tf) {
    if (!(tf %in% rownames(dom@features))) {
        stop("Transcription factor not found in domino object")
    }
    test_dom <- dom
    test_dom@features[tf, ] <- 0
    return(test_dom)
}