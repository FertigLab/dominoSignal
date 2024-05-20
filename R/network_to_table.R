# For additional helper functions that might prove useful in any domino2 analysis
# (with considerable assistance from Sushma and others from the Fertig lab)

require(purrr)
require(tidyverse)
require(Matrix)

#' Convert between ligand names and gene names
#' @param dom A built domino object
#' @param genes A vector of genes on which to resolve ligand and gene names
#' @return A vector of names where ligand names have been replaced with gene names if applicable
#' @keywords internal
resolve_names <- function(dom, genes) {
    # Grab rl_map from inside object
    rl_map <- dom@misc[["rl_map"]]
    # For each gene supplied
    genes_resolved <- sapply(genes, function(l) {
        # For each ligand, find the first row in the rl_map
        int <- rl_map[rl_map$L.name == l, ][1, ]
        # If the ligand is in the rl_map and ligand name and ligand gene aren't the same,
        # and there is no comma in the ligand gene, use the gene
        if ((l %in% rl_map$L.name) & (int$L.name != int$L.gene) & !grepl("\\,", int$L.gene)) {
            remove_semicolons(int$L.gene)
            # Otherwise, use the ligand name
        } else if (l %in% rl_map$L.name) {
            remove_semicolons(int$L.name)
            # If the ligand is not in the rl_map, just return ligand without semicolons
        } else {
            remove_semicolons(l)
        }
    })
    return(genes_resolved)
}

#' Convert between complex names and gene names
#' @param dom A built domino object
#' @param genes A vector of genes on which to resolve ligand and gene names
#' @return A list of names where ligand names have been replaced with gene names if applicable
#' @keywords internal
resolve_complexes <- function(dom, genes) {
    # For each gene, make a list such that
    genes_list <- lapply(genes, function(l) {
        # If ligand is in a complex
        if (length(l) > 1) {
            # Just check for semis and return? Why is this happening?
            check_semis <- sapply(X = l, FUN = remove_semicolons, USE.NAMES = FALSE)
        } else if (l %in% names(dom@linkages$complexes)) {
            # Return the genes in the complex (after checking for semicolons)
            subunits = dom@linkages$complexes[[l]]
            check_semis <- sapply(X = subunits, FUN = remove_semicolons, USE.NAMES = FALSE)
            # Otherwise just return the ligand (after checking for semicolons)
        } else {
            check_semis <- remove_semicolons(l)
        }
    })
    # Name the list of gene(s) with the original names
    names(genes_list) <- names(genes)
    return(genes_list)
}

#' Get average expression for complexes
#' @param exp_mat A matrix of genes x clusters, where the values are
#'  z_scores averaged over the clusters
#' @param complexes_list a list like the one at dom@linkages$complexes
#' @return A list of average expression for each complex
#' @keywords internal
avg_exp_for_complexes <- function(exp_mat, complexes_list) {
    # Trim the complexes list to only include those with genes in the data
    trim_list <- complexes_list %>% keep(~ {
        all(.x %in% rownames(exp_mat))
    })
    # For each included item, create a list
    gene_exp_list <- lapply(seq_along(trim_list), function(x) {
        # If there is more than 1 item in the complex
        if (length(trim_list[[x]]) > 1) {
            # Return the average expression of those items
            return(colMeans(exp_mat[unlist(trim_list[[x]]), ]))
        } else {
            # If there's 1 item, no need to average
            return(exp_mat[unlist(trim_list[[x]]), ])
        }
    })
    # Name the results based on the included complexes
    names(gene_exp_list) <- names(trim_list)
    return(gene_exp_list)
}

#' Get ligands with resolved names
#' @param dom A built domino object
#' @return A list of ligands and complexes with resolved names
#' @keywords internal
get_resolved_ligands <- function(dom) {
    all_lig <- unlist(dom@linkages$rec_lig)
    all_lig <- unique(all_lig)
    all_lig <- all_lig[!all_lig == ""]
    all_lig_names_resolved <- resolve_names(dom, all_lig)
    
    if (length(dom@linkages$complexes) > 0) {
        all_lig_complexes_resolved_list <- resolve_complexes(dom, all_lig_names_resolved)
        all_lig_names_resolved <- unlist(all_lig_complexes_resolved_list)
    }
    
    return(list("lig_names" = unique(all_lig_names_resolved), "complex_names" = all_lig_complexes_resolved_list))
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
    cl_ligands <- matrix(0, ncol = length(send_clusters), nrow = length(lig_genes))
    colnames(cl_ligands) <- send_clusters
    rownames(cl_ligands) <- lig_genes
    
    for (c2 in send_clusters) {
        n_cell <- length(which(dom@clusters == c2))
        if (n_cell > 1) {
            if (exp_type == "counts") {
                sig <- rowMeans(dom@counts[lig_genes, which(dom@clusters == c2)])
            } else if (exp_type == "z_scores") {
                sig <- rowMeans(dom@z_scores[lig_genes, which(dom@clusters == c2)])
            }
        } else if (n_cell == 1) {
            if (exp_type == "counts") {
                sig <- dom@counts[lig_genes, which(dom@clusters == c2)]
            } else if (exp_type == "z_scores") {
                sig <- dom@z_scores[lig_genes, which(dom@clusters == c2)]
            }
            
        } else {
            sig <- rep(NA, length(lig_genes))
            names(sig) <- lig_genes
        }
        
        cl_ligands[, c2] <- sig
    }    
    if (length(dom@linkages$complexes) > 0) {
        cl_ligands_coll_list <- avg_exp_for_complexes(cl_ligands, complexes)
        
        if (length(cl_ligands_coll_list) > 1) {
            cl_ligands <- do.call(rbind, cl_ligands_coll_list)
        }
    }
    
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
    row_list <- list()
    i = 1
    for (cl in rec_clusters) {
        for (tf in names(dom@linkages$clust_tf_rec[[cl]])) {
            recs <- dom@linkages$clust_tf_rec[[cl]][[tf]]
            
            if (length(recs)) {
                for (rec in recs) {
                    rec_cell <- length(which(dom@clusters == cl))
                    rec_sep = unlist(resolve_complexes(dom, rec))
                    
                    if (rec_cell > 1) {
                        if (length(rec_sep) > 1) {
                            if (exp_type == "counts") {
                                rec_sig <- mean(rowMeans(dom@counts[rec_sep, which(dom@clusters == cl)]))
                            } else if (exp_type == "z_scores") {
                                rec_sig <- mean(rowMeans(dom@z_scores[rec_sep, which(dom@clusters == cl)]))
                            }
                        } else {
                            if (exp_type == "counts") {
                                rec_sig <- mean(dom@counts[rec_sep, which(dom@clusters == cl)])
                            } else if (exp_type == "z_scores") {
                                rec_sig <- mean(dom@z_scores[rec_sep, which(dom@clusters == cl)])
                            }
                        }
                        tf_sig <- mean(dom@features[tf, which(dom@clusters == cl)])
                    } else if (rec_cell == 1) {
                        if (length(rec_sep) > 1) {
                            if (exp_type == "counts") {
                                rec_sig <- mean(dom@counts[rec_sep, which(dom@clusters == cl)])
                            } else if (exp_type == "z_scores") {
                                rec_sig <- mean(dom@z_scores[rec_sep, which(dom@clusters == cl)])
                            }
                        } else {
                            if (exp_type == "counts") {
                                rec_sig <- dom@counts[rec_sep, which(dom@clusters == cl)]
                            } else if (exp_type == "z_scores") {
                                rec_sig <- dom@z_scores[rec_sep, which(dom@clusters == cl)]
                            }
                        }
                        tf_sig <- dom@features[tf, which(dom@clusters == cl)]
                    } else {
                        rec_sig <- NA
                        tf_sig <- NA
                    }
                    
                    ligs <- dom@linkages$rec_lig[[rec]]
                    ligs <- resolve_names(dom, ligs)
                    
                    if (length(ligs)) {
                        df_temp <- cl_ligands_sub[cl_ligands_sub$ligand %in% ligs, ]
                        
                        if (nrow(df_temp)) {
                            new_row <- data.frame(
                                ligand = df_temp$ligand,
                                receptor = rec,
                                transcription_factor = tf,
                                ligand_exp = df_temp$mean_counts,
                                rec_exp = rec_sig,
                                tf_auc = tf_sig,
                                sending_cl = df_temp$cluster,
                                receiving_cl = cl
                            )
                            row_list[[i]] <- new_row
                            i = i + 1
                        }
                    }
                }
            }
        }
    }
    df = purrr::list_rbind(row_list)
    return(df)
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
#' data("pbmc_dom_tiny_built")
#' df <- network_to_df(pbmc_dom_tiny_built, exp_type = "z_scores")
#' head(df)
network_to_df <- function(dom, send_clusters = NULL, rec_clusters = NULL, exp_type = c("counts", "z_scores")) {
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
    
    df <- get_signaling_info(dom, rec_clusters, cl_ligands_sub, exp_type)
    
    return(df)
}

#' A helper function for removing semicolons
#' @param gene A gene to be checked for a semicolon (and removed if present)
#' @return The gene with semicolons removed
#' @keywords internal
remove_semicolons <- function(gene) {
    semi_removed = gsub(";^", "", gene, fixed = TRUE)
    if (stringr::str_count(semi_removed, pattern = ";")) {
        sep_genes = unlist(strsplit(semi_removed, ";"))
    } else {
        sep_genes = semi_removed
    }
    
    return(sep_genes)
}
