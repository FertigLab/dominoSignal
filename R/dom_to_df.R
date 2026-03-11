#' @importFrom reshape2 melt
#' @importFrom purrr list_rbind
NULL


#' Get ligands with resolved names
#' @param dom A built domino object
#' @return A list of ligands and complexes with resolved names
#' @keywords internal
get_resolved_ligands <- function(dom) {
    
    check_arg(dom, allow_class = "domino", allow_len = 1)

    all_lig <- unlist(dom@linkages$rec_lig)
    all_lig <- unique(all_lig)
    all_lig <- all_lig[nzchar(all_lig, keepNA = TRUE)]
    all_lig_names_resolved <- resolve_names(dom, all_lig)

    if (length(dom@linkages$complexes) > 0) {
        lig_complexes_resolved_list <- resolve_complexes(dom, all_lig_names_resolved)
        all_lig_names_resolved <- unlist(lig_complexes_resolved_list)
    }

    return(list("lig_names" = unique(all_lig_names_resolved), "complex_names" = lig_complexes_resolved_list))
}

#' Get ligand expression matrix for outgoing clusters
#' @param dom A built domino object
#' @param send_clusters A cluster name or vector for which the outgoing signal is desired
#' @param lig_genes A vector of ligand genes
#' @param complexes A list of complexes with names as complex and values as component genes
#' @param exp_type A character vector of length 1, either "counts" or "z_scores"
#' @return A matrix of ligand expression for outgoing clusters
#' @keywords internal
get_ligand_expression <- function(dom, send_clusters, lig_genes, complexes, exp_type) {
    
    check_arg(dom, allow_class = "domino", allow_len = 1)
    check_arg(send_clusters, allow_class = c("character", "factor"), allow_values = dom_clusters(dom))
    check_arg(lig_genes, allow_class = "character")
    check_arg(complexes, allow_class = "list", need_names = TRUE)
    check_arg(exp_type, allow_class = "character", allow_values = c("counts", "z_scores"), allow_len = 1)

    expr_mat <- switch(exp_type,
        counts = dom_counts(dom),
        z_scores = dom_zscores(dom),
        stop("Invalid exp_type. Must be either 'counts' or 'z_scores'.")
    )

    cl_ligands <- matrix(NA_real_, nrow = length(lig_genes), ncol = length(send_clusters),
        dimnames = list(lig_genes, send_clusters))

    for (clust in send_clusters) {
        idx <- which(dom@clusters == clust)
        if (length(idx) == 0) next
        cl_ligands[ , clust] <- rowMeans(expr_mat[lig_genes, idx, drop = FALSE])
    }

    if (length(dom_linkages(dom, "complexes")) > 0) {
        cl_ligands_coll_list <- avg_exp_for_complexes(cl_ligands, complexes)
        if (length(cl_ligands_coll_list) > 0) {
            cl_ligands <- Reduce(rbind, cl_ligands_coll_list)
        }
    }

    cl_ligands <- as.matrix(cl_ligands)

    return(cl_ligands)
}

#' Get ligand-receptor signaling information
#' @param dom A built domino object
#' @param rec_clusters A cluster name or vector for which incoming signal is desired
#' @param cl_ligands_sub A data frame of ligand expression for outgoing clusters
#' @param exp_type A character vector of length 1, either "counts" or "z_scores"
#' @return A data frame of signaling information
#' @keywords internal
get_signaling_info <- function(dom, rec_clusters, cl_ligands_sub, exp_type) {

    check_arg(dom, allow_class = "domino", allow_len = 1)
    check_arg(rec_clusters, allow_class = c("character", "factor"), allow_values = dom_clusters(dom))
    check_arg(cl_ligands_sub, allow_class = "data.frame")
    check_arg(exp_type, allow_class = "character", allow_values = c("counts", "z_scores"), allow_len = 1)

    expr_mat <- switch(exp_type,
        counts = dom_counts(dom),
        z_scores = dom_zscores(dom),
        stop("Invalid exp_type. Must be either 'counts' or 'z_scores'.")
    )

    tf_mat <- dom_tf_activation(dom)
    row_list <- list()

    for (cl in rec_clusters) {
        for (tf in names(dom@linkages$clust_tf_rec[[cl]])) {
            recs <- dom@linkages$clust_tf_rec[[cl]][[tf]]

            if (length(recs) == 0) next

            rec_idx <- which(dom@clusters == cl)

            if (length(rec_idx) > 0) {
                tf_sig <- mean(tf_mat[tf, rec_idx, drop = FALSE])
            } else {
                tf_sig <- NA
            }

            for (rec in recs) {
                rec_sep <- unlist(resolve_complexes(dom, rec))

                if (length(rec_idx) > 0) {
                    rec_sig <- mean(expr_mat[rec_sep, rec_idx, drop = FALSE])
                } else {
                    rec_sig <- NA
                }

                ligs <- resolve_names(dom, dom@linkages$rec_lig[[rec]])

                if (length(ligs) > 0) {
                    df_tmp <- cl_ligands_sub[cl_ligands_sub$ligand %in% ligs, ]

                    if (nrow(df_tmp) > 0) {
                        new_row <- data.frame(
                            ligand = df_tmp$ligand,
                            receptor = rec,
                            transcription_factor = tf,
                            ligand_exp = df_tmp$mean_counts,
                            rec_exp = rec_sig,
                            tf_auc = tf_sig,
                            sending_cl = df_tmp$cluster,
                            receiving_cl = cl
                        )
                        row_list <- c(row_list, list(new_row))
                    }
                }
            }
        }
    }

    if (length(row_list) == 0) {
        return(data.frame())
    }       

    dframe <- purrr::list_rbind(row_list)
    return(dframe)
}

#' Turn domino object signaling information into a data frame
#'
#' This function takes a domino object and returns a data frame of signaling information with
#' columns for ligand, receptor, transcription factor, ligand expression, receptor expression,
#' transcription factor expression, sending cluster, and receiving cluster. Expression values are
#' averaged counts across cells in a given cluster (sending for the ligand, receiving for the receptor
#' and transcription factor).
#'
#' @param dom A domino object
#' @param send_clusters Cluster name(s) for which the outgoing signal is desired
#' @param rec_clusters Cluster name(s) for which incoming signal is desired
#' @param exp_type A character vector of length 1, either "counts" or "z_scores", to indicate desired
#'  expression type
#' @return A data.frame of signaling information, with columns for ligand, receptor,
#'  transcription factor, ligand expression, receptor expression,
#'  transcription factor expression, sending cluster, and receiving cluster
#' @export
#' @examples
#' data("DominoObjects")
#' dom <- DominoObjects$built_dom_tiny
#' df <- dom_to_df(dom, exp_type = "z_scores")
#' head(df)
dom_to_df <- function(dom, send_clusters = NULL, rec_clusters = NULL, exp_type = c("counts", "z_scores")) {

    check_arg(dom, allow_class = "domino", allow_len = 1)
    if (!is.null(send_clusters)) {
        check_arg(send_clusters, allow_class = "character", allow_values = dom_clusters(dom))
    }
    if (!is.null(rec_clusters)) {
        check_arg(rec_clusters, allow_class = "character", allow_values = dom_clusters(dom))
    }
    check_arg(exp_type, allow_class = "character", allow_values = c("counts", "z_scores"), allow_len = 1)

    all_lig_names_resolved_lists <- get_resolved_ligands(dom)
    all_lig_names_resolved <- all_lig_names_resolved_lists$lig_names
    all_lig_complexes_resolved <- all_lig_names_resolved_lists$complex_names

    if (exp_type == "counts") {
        lig_genes <- intersect(all_lig_names_resolved, rownames(dom@counts))
    } else if (exp_type == "z_scores") {
        lig_genes <- intersect(all_lig_names_resolved, rownames(dom@z_scores))
    }

    if (is.null(send_clusters)) {
        send_clusters <- levels(dom@clusters)
    }

    cl_ligands <- get_ligand_expression(dom, send_clusters, lig_genes, all_lig_complexes_resolved, exp_type)

    cl_ligands_sub <- reshape2::melt(cl_ligands)
    colnames(cl_ligands_sub) <- c("ligand", "cluster", "mean_counts")

    if (is.null(rec_clusters)) {
        rec_clusters <- levels(dom@clusters)
    }

    dframe <- get_signaling_info(dom, rec_clusters, cl_ligands_sub, exp_type)

    return(dframe)
}
