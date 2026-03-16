#' @importFrom reshape2 melt
#' @importFrom purrr list_rbind
NULL


#' Get ligands with resolved names
#' @param dom A built domino object (as output by [build_domino()])
#' @return A named list with two elements:
#'   \item{lig_names}{Character vector of unique ligand gene names or aliases}
#'   \item{complex_names}{Named list mapping complex names to component gene vectors}
#' @keywords internal
get_resolved_ligands <- function(dom) {
    
    check_arg(dom, allow_class = "domino", allow_len = 1)

    all_lig <- unlist(dom@linkages$rec_lig)
    all_lig <- unique(all_lig)
    all_lig <- all_lig[nzchar(all_lig)]
    all_lig_names_resolved <- resolve_names(dom, all_lig)

    lig_complexes_resolved_list <- resolve_complexes(dom, all_lig_names_resolved)
    all_lig_names_resolved <- unlist(lig_complexes_resolved_list)

    return(list("lig_names" = unique(all_lig_names_resolved), "complex_names" = lig_complexes_resolved_list))
}

#' Get ligand expression matrix for outgoing clusters
#' @param dom A built domino object (as output by [build_domino()])
#' @param send_clusters Character or factor vector of cluster names for which to compute outgoing ligand signals
#' @param lig_genes Character vector of ligand gene names; should intersect with expression matrix rownames
#' @param complexes Named list where names are complex identifiers and values are character vectors of component genes.
#'   Used to aggregate expression for multi-gene complexes. Empty list is acceptable.
#'   The output of [get_resolved_ligands()] can be used here.
#' @param exp_type Character of length 1: either "counts" or "z_scores" to specify expression type
#' @return Matrix with ligands/complexes as rows and send_clusters as columns; entries are mean expression.
#'   Missing ligands and empty clusters yield NA; row and column names are preserved.
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

    if (!all(lig_genes %in% rownames(expr_mat))) {
        message("Some ligands not found in expression matrix: ",
            toString(setdiff(lig_genes, rownames(expr_mat))))
        lig_genes <- intersect(lig_genes, rownames(expr_mat))
        if (length(lig_genes) == 0) {
            stop("No ligands found in expression matrix")
        }
    }

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
            rownames(cl_ligands) <- names(cl_ligands_coll_list)
        }
    }

    cl_ligands <- as.matrix(cl_ligands)

    return(cl_ligands)
}

#' Get ligand-receptor signaling information
#' @param dom A built domino object (as output by [build_domino()])
#' @param rec_clusters Character or factor vector of cluster names for which to compute incoming receptor/TF signals
#' @param cl_ligands_sub Data frame with columns 'ligand', 'cluster', 'mean_counts';
#'  using [reshape2::melt()] on the output of [get_ligand_expression()] is a convenient way to get this
#' @param exp_type Character of length 1: either "counts" or "z_scores" to specify expression type
#' @return Data frame with columns:
#'   ligand, receptor, transcription_factor, ligand_exp, rec_exp, tf_auc, sending_cl, receiving_cl.
#'   Each row represents a ligand-receptor-TF triplet with receiver cluster, and corresponding expression values.
#'   Returns empty data frame if no interactions found.
#' @keywords internal
get_signaling_info <- function(dom, rec_clusters, cl_ligands_sub, exp_type) {

    check_arg(dom, allow_class = "domino", allow_len = 1)
    check_arg(rec_clusters, allow_class = c("character", "factor"), allow_values = dom_clusters(dom))
    check_arg(cl_ligands_sub, allow_class = "data.frame", need_vars = c("ligand", "cluster", "mean_counts"))
    check_arg(exp_type, allow_class = "character", allow_values = c("counts", "z_scores"), allow_len = 1)

    expr_mat <- switch(exp_type,
        counts = dom_counts(dom),
        z_scores = dom_zscores(dom),
        stop("Invalid exp_type. Must be either 'counts' or 'z_scores'.")
    )

    tf_mat <- dom_tf_activation(dom)
    row_list <- list()

    # Helper function to return NA for empty indices and mean where available
    mean_or_na <- function(mat, rows, cols) {
        if (length(cols) == 0) return(NA_real_)
        mean(mat[rows, cols, drop = FALSE])
    }

    for (cl in rec_clusters) {
        rec_idx <- which(dom@clusters == cl)
        tf_idx <- dom@linkages$clust_tf_rec[[cl]]
        if (length(tf_idx) == 0) next
        for (tf in names(tf_idx)) {
            recs <- tf_idx[[tf]]
            if (length(recs) == 0) next

            tf_sig <- mean_or_na(tf_mat, tf, rec_idx)

            for (rec in recs) {
                rec_sep <- unlist(resolve_complexes(dom, rec))
                rec_sig <- mean_or_na(expr_mat, rec_sep, rec_idx)

                ligs <- resolve_names(dom, dom@linkages$rec_lig[[rec]])
                if (length(ligs) == 0) next

                df_tmp <- cl_ligands_sub[cl_ligands_sub$ligand %in% ligs, ]
                if (nrow(df_tmp) == 0) next

                row_list[[length(row_list) + 1]] <- data.frame(
                    ligand = df_tmp$ligand,
                    receptor = rec,
                    transcription_factor = tf,
                    ligand_exp = df_tmp$mean_counts,
                    rec_exp = rec_sig,
                    tf_auc = tf_sig,
                    sending_cl = df_tmp$cluster,
                    receiving_cl = cl
                )
            }
        }
    }

    if (length(row_list) == 0) {
        message("No interactions found for the specified clusters and expression type.")
        return(data.frame())
    }       

    return(purrr::list_rbind(row_list))
}

#' Turn domino object signaling information into a data frame
#'
#' Constructs a data frame of ligand-receptor-TF signaling triplets from a built domino object.
#' Expression values are averaged across cells within each cluster (ligands from sending clusters,
#' receptors and TFs from receiving clusters).
#'
#' @param dom A built domino object (as output by [build_domino()])
#' @param send_clusters Character/factor vector of cluster names for ligand signals.
#'   If NULL (default), uses all clusters in dom.
#' @param rec_clusters Character/factor vector of cluster names for receptor/TF signals.
#'   If NULL (default), uses all clusters in dom.
#' @param exp_type Character of length 1: either "counts" or "z_scores" for expression type
#' @return Data frame with columns: ligand, receptor, transcription_factor,
#'   ligand_exp, rec_exp, tf_auc, sending_cl, receiving_cl.
#'   Each row is a ligand-receptor-TF triplet from a sender-to-receiver cluster pair with
#'   corresponding mean expression values.
#'   Empty data frame returned if no valid interactions found.
#' @export
#' @examples
#' data("DominoObjects")
#' dom <- DominoObjects$built_dom_tiny
#' df <- dom_to_df(dom, exp_type = "z_scores")
#' head(df)
dom_to_df <- function(dom, send_clusters = NULL, rec_clusters = NULL, exp_type = c("counts", "z_scores")) {

    check_arg(dom, allow_class = "domino", allow_len = 1)
    if (!is.null(send_clusters)) {
        check_arg(send_clusters, allow_class = c("character", "factor"), allow_values = dom_clusters(dom))
    }
    if (!is.null(rec_clusters)) {
        check_arg(rec_clusters, allow_class = c("character", "factor"), allow_values = dom_clusters(dom))
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
