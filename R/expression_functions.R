#' Get average expression for complexes
#'
#' @param exp_mat A matrix(or dataframe) of genes x clusters, values are z-scores averaged over the clusters
#' @param complexes_list A list similar to dom@linkages$complexes
#'
#' @return A list containing average expression for any complexes
#' @keywords internal
avg_exp_for_complexes <- function(exp_mat, complexes_list) {
    check_arg(exp_mat, allow_class = c("matrix", "data.frame"), need_rownames = TRUE, need_colnames = TRUE)
    check_arg(complexes_list, allow_class = "list")
    # dplyr verbs require a data frame; coerce here so callers can pass either
    exp_mat <- as.data.frame(exp_mat)
    # Trim the complexes list to only include those with genes in the data
    trim_list <- complexes_list |>
        purrr::keep(~ {
            all(.x %in% rownames(exp_mat))
        })
    gene_exp_list <- lapply(seq_along(trim_list), function(x) {
        if (length(trim_list[[x]]) > 1) {
            mean_exp <- exp_mat |>
                dplyr::filter(rownames(exp_mat) %in% trim_list[[x]]) |>
                dplyr::summarise(dplyr::across(dplyr::all_of(colnames(exp_mat)), mean))
            return(mean_exp)
        } else {
            return(exp_mat[trim_list[[x]], , drop = FALSE])
        }
    })
    names(gene_exp_list) <- names(trim_list)
    return(gene_exp_list)
}

#' Get average expression for a set of genes over cluster(s)
#'
#' @param dom A domino object
#' @param clusts Cluster(s) for which we want to get average expression
#' @param genes The genes for which we want to get average expression
#'
#' @return A dataframe of genes x clusters, values are z-scores averaged over the clusters
#' @keywords internal
mean_exp_by_cluster <- function(dom, clusts, genes) {
    check_arg(dom, allow_class = "domino", allow_len = 1)
    check_arg(clusts, allow_class = "character")
    check_arg(genes, allow_class = "character")

    gene_exp_list <- purrr::map(seq_along(clusts), function(x) {
        cl <- clusts[x]
        n_cell <- length(which(dom@clusters == cl))
        if (n_cell > 1) {
            sig <- rowMeans(dom@z_scores[genes, which(dom@clusters == cl)])
        } else if (n_cell == 1) {
            sig <- dom@z_scores[genes, which(dom@clusters == cl)]
        } else {
            sig <- rep(0, length(genes))
            names(sig) <- genes
        }
        sig[which(sig < 0)] <- 0
        sig <- as.data.frame(sig)
        colnames(sig) <- cl
        return(sig)
    })
    gene_exp <- purrr::list_cbind(gene_exp_list)
    return(gene_exp)
}

#' Calculate mean ligand expression as a data frame for plotting in circos plot
#'
#' Creates a data frame of mean ligand expression for use in plotting a circos
#' plot of ligand expression and saving tables of mean expression.
#'
#' @param x Gene by cell expression matrix
#' @param ligands Character vector of ligand genes to be quantified
#' @param cell_ident Vector of cell type (identity) names for which to calculate mean ligand gene expression
#' @param cell_barcodes Vector of cell barcodes (colnames of x) belonging to cell_ident to calculate mean
#'   expression across
#' @param destination Name of the receptor with which each ligand interacts
#' @return A data frame of ligand expression targeting the specified receptor
#' @export
#' @examples
#' data(DominoObjects)
#' counts <- dom_counts(DominoObjects$built_dom_tiny)
#' mean_exp <- mean_ligand_expression(counts,
#'     ligands = c("PTPRC", "FASLG"), cell_ident = "CD14_monocyte",
#'     cell_barcodes = colnames(counts), destination = "FAS"
#' )
#'
mean_ligand_expression <- function(x, ligands, cell_ident, cell_barcodes, destination) {
    check_arg(x, allow_class = c("matrix", "data.frame"))
    check_arg(ligands, allow_class = "character")
    check_arg(cell_ident, allow_class = "character")
    check_arg(cell_barcodes, allow_class = "character")
    check_arg(destination, allow_class = "character", allow_len = 1)

    # initiate data frame to store results
    dframe <- NULL
    for (feat in ligands) {
        # index of ligand row
        lig_index <- grep(paste0("^", feat, "$"), rownames(x))
        # column indecies of cells belonging to cell_ident
        cell_index <- colnames(x) %in% cell_barcodes
        cell_df <- data.frame(
            origin = paste0(cell_ident, "_", feat),
            destination = destination,
            mean.expression = mean(x[lig_index, cell_index])
        )
        dframe <- rbind(dframe, cell_df)
    }
    return(dframe)
}

#' Normalize a matrix to its max value by row or column
#'
#' Normalizes a matrix to its max value by row or column
#'
#' @param mat Matrix to be normalized
#' @param dir Direction to normalize the matrix (either "row" for row or "col" for column)
#' @return A normalized matrix in the direction specified.
#' @keywords internal
#'
do_norm <- function(mat, dir) {
    check_arg(mat, allow_class = c("matrix", "data.frame"))
    check_arg(dir, allow_class = "character", allow_len = 1, allow_values = c("row", "col"))

    if (dir == "row") {
        mat <- t(apply(mat, 1, function(x) {
            x / max(x)
        }))
        return(mat)
    } else if (dir == "col") {
        mat <- apply(mat, 2, function(x) {
            x / max(x)
        })
        return(mat)
    }
}
