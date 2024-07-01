# A merge function for comparing/combining objects

#' Merge db_info slot of domino objects
#' 
#' @param dom_list A list of domino objects to be merged
#' @keywords internal
merge_db_info <- function(dom_list) { # This is going to break with non-cellphonedb inputs
    db_info_list = lapply(dom_list, slot, "db_info")
    # Get union of matrices in db_info_list using dplyr::union
    db_info_merged = reduce(db_info_list, dplyr::union)
    return(db_info_merged)
}

#' Merge build and create items in misc slot
#' 
#' @param misc_list A list of misc items in domino objects to be merged
#' @keywords internal
merge_misc_build_create <- function(misc_list) {
    # Get list item named build from each item in misc_list
    build_list = lapply(misc_list, function(x) {
        x$build
    })
    builds = unique(unlist(build_list))
    create_list = lapply(misc_list, function(x) {
        x$create
    })
    creates = unique(unlist(create_list))
    return(list("build" = builds, "create" = creates))
}

#' Merge rl_map slot of misc slot of domino objects
#' 
#' @param misc_list A list of misc items in domino objects to be merged
#' @keywords internal
merge_misc_rl_map <- function(misc_list) {
    rl_map_list = lapply(misc_list, function(x) {
        x["rl_map"]$rl_map
    })
    # Get union of matrices in rl_map_list using dplyr::union
    rl_map_merged = reduce(rl_map_list, dplyr::union)
    return(rl_map_merged)
}

#' Merge cor slot of misc slot of domino objects
#' 
#' @param misc_list A list of misc items in domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_misc_cor <- function(misc_list, value_method) {
    cor_list = lapply(misc_list, function(x) {
        x["rec_cor"]$rec_cor
    })
    # Get union of rownames and column names in cor_list
    row_cor = unique(unlist(lapply(cor_list, rownames)))
    col_cor = unique(unlist(lapply(cor_list, colnames)))
    # Get mean value of each combination of row and column names
    cor_merged = matrix(0, nrow = length(row_cor), ncol = length(col_cor))
    rownames(cor_merged) = row_cor
    colnames(cor_merged) = col_cor
    for (r in row_cor) {
        for (c in col_cor) {
            if (value_method == "max") {
                cor_merged[r, c] = max(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "min") {
                cor_merged[r, c] = min(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "mean") {
                cor_merged[r, c] = mean(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "range") {
                cor_merged[r, c] = max(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                }))) - min(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else {
                stop("Invalid method for value merging")
            }
        }
        return(cor_merged)
    }

}

#' Merge cl_rec_percent slot of misc slot of domino objects
#' 
#' @param misc_list A list of misc items in domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_misc_percent <- function(misc_list, value_method) {
    percent_list = lapply(misc_list, function(x) {
        x["cl_rec_percent"]$cl_rec_percent
    })
    # Same as merge_misc_cor but for percent
    row_percent = unique(unlist(lapply(percent_list, rownames)))
    col_percent = unique(unlist(lapply(percent_list, colnames))) # Is this problematic if there are different clusters?
    percent_merged = matrix(0, nrow = length(row_percent), ncol = length(col_percent))
    rownames(percent_merged) = row_percent
    colnames(percent_merged) = col_percent
    for (r in rownames(percent_merged)) {
        for (c in colnames(percent_merged)) {
            if (value_method == "max") {
                percent_merged[r, c] = max(unlist(lapply(percent_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "min") {
                percent_merged[r, c] = min(unlist(lapply(percent_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "mean") {
                percent_merged[r, c] = mean(unlist(lapply(percent_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "range") {
                percent_merged[r, c] = max(unlist(lapply(percent_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE) - min(unlist(lapply(percent_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)   
            } else {
                stop("Invalid method for value merging")
            }
        }
    }
    return(percent_merged)
}

#' Merge build_vars slot of misc slot of domino objects
#' 
#' @param misc_list A list of misc items in domino objects to be merged
#' @keywords internal
merge_misc_build <- function(misc_list) { # May need to change structure in object definition to list
    build_vars = lapply(misc_list, function(x) {
        x["build_vars"]$build_vars
    })
    # For each variable in build_vars, see if there is more than one unique value
    merged = list()
    for (var in seq_along(names(build_vars[[1]]))) {
        if (length(unique(unlist(lapply(build_vars, function(x) {x[[var]]})))) > 1) {
            # Get that variable from each item in build_vars
            var_list = lapply(build_vars, function(x) {
                x[[var]]
            })
            # Add object name to each item in var_list
            merged[[names(build_vars[[1]][var])]] = unlist(var_list)
        } else {
            # If there is only one unique value, just use that
            merged[[names(build_vars[[1]][var])]] = build_vars[[1]][[var]]
        }
    }
    return(merged)
}

#' Merge misc slot of domino objects
#' 
#' Merge the lists of misc items in domino objects
#' 
#' @param dom_list A list of domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_misc <- function(dom_list, value_method) {
    misc_list = lapply(dom_list, slot, "misc")
    merged_misc = list(
        "build" = merge_misc_build_create(misc_list)[["build"]],
        "create" = merge_misc_build_create(misc_list)[["create"]],
        "rl_map" = merge_misc_rl_map(misc_list),
        "rec_cor" = merge_misc_cor(misc_list, value_method),
        "cl_rec_percent" = merge_misc_percent(misc_list, value_method),
        "build_vars" = merge_misc_build(misc_list)
    )
    return(merged_misc)
}

#' Merge counts slot of domino objects
#' 
#' Merge list of counts matrices and select features if desired
#' 
#' @param dom_list A list of domino objects to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_counts <- function(dom_list, feature_method) {
    counts_list = lapply(dom_list, slot, "counts")
    # Check whether there are any duplicate column names
    if (length(unique(unlist(lapply(counts_list, colnames)))) != length(unlist(lapply(counts_list, colnames)))) {
        # If not, add name of domino object to each column name as prefix
        counts_list = lapply(seq_along(counts_list), function(x) {
            colnames(counts_list[[x]]) = paste0(names(dom_list)[x], "_", colnames(counts_list[[x]]))
            return(counts_list[[x]])
        })
    }
    # Once we have unique column names, we can merge the matrices, using RVenn to select features
    all_feats = lapply(counts_list, rownames)
    all_feats = RVenn::Venn(all_feats)
    if (feature_method == "intersect") {
        feats = RVenn::overlap(all_feats)
    } else if (feature_method == "union") {
        feats = RVenn::unite(all_feats)
    } else if (feature_method == "outersect") {
        # Use RVenn::discern to find features unique to each item in all_feats
        feats = c()
        for (i in seq_along(all_feats)) {
            unique_items = RVenn::discern(all_feats, i)
            feats = c(feats, unique_items)
        }
    } else {
        stop("Invalid method for feature merging")
    }
    merged = do.call(cbind, counts_list)
    merged = merged[feats, ]
    return(merged)
}

#' Merge clust_de slot from domino objects
#' 
#' Using method, merge DE results from domino objects
#' 
#' @param dom_list A list of domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_de <- function(dom_list, value_method) {
    de_list = lapply(dom_list, slot, "clust_de")

    row_de = unique(unlist(lapply(de_list, rownames)))
    col_de = unique(unlist(lapply(de_list, colnames))) # Is this problematic if there are different clusters?
    de_merged = matrix(0, nrow = length(row_de), ncol = length(col_de))
    rownames(de_merged) = row_de
    colnames(de_merged) = col_de
    for (r in rownames(de_merged)) {
        for (c in colnames(de_merged)) {
            if (value_method == "max") {
                de_merged[r, c] = max(unlist(lapply(de_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "min") {
                de_merged[r, c] = min(unlist(lapply(de_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "mean") {
                de_merged[r, c] = mean(unlist(lapply(de_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else if (value_method == "range") {
                de_merged[r, c] = max(unlist(lapply(de_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE) - min(unlist(lapply(de_list, function(x) {
                    rc_val_mat(x, r, c)
                })), na.rm = TRUE)
            } else {
                stop("Invalid method for value merging")
            }
        }
    }
    return(de_merged)
}

#' Merge z_scores slot of domino objects
#' 
#' Merge list of z_scores matrices and select features if desired
#' 
#' @param dom_list A list of domino objects to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_zs <- function(dom_list, feature_method) {
    zs_list = lapply(dom_list, slot, "z_scores")
    # Check whether there are any duplicate column names
    if (length(unique(unlist(lapply(zs_list, colnames)))) != length(unlist(lapply(zs_list, colnames)))) {
        # If not, add name of domino object to each column name as prefix
        zs_list = lapply(seq_along(zs_list), function(x) {
            colnames(zs_list[[x]]) = paste0(names(dom_list)[x], "_", colnames(zs_list[[x]]))
            return(zs_list[[x]])
        })
    }
    # Once we have unique column names, we can merge the matrices, using RVenn to select features
    all_feats = lapply(zs_list, rownames)
    all_feats = RVenn::Venn(all_feats)
    if (feature_method == "intersect") {
        feats = RVenn::overlap(all_feats)
    } else if (feature_method == "union") {
        feats = RVenn::unite(all_feats)
    } else if (feature_method == "outersect") {
        feats = c()
        for (i in seq_along(zs_list)) {
            unique_items = RVenn::discern(all_feats, i)
            feats = c(feats, unique_items)
        }
    } else {
        stop("Invalid method for feature merging")
    }
    merged = do.call(cbind, zs_list)
    merged = merged[feats, ]
    return(merged)

}

#' Merge correlation slot of domino objects
#' 
#' Merge list of correlation matrices
#' 
#' @param dom_list A list of domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_cor <- function(dom_list, value_method) {
    cor_list = lapply(dom_list, slot, "cor")
    row_cor = unique(unlist(lapply(cor_list, rownames)))
    col_cor = unique(unlist(lapply(cor_list, colnames)))
    cor_merged = matrix(0, nrow = length(row_cor), ncol = length(col_cor))
    rownames(cor_merged) = row_cor
    colnames(cor_merged) = col_cor
    for (r in row_cor) {
        for (c in col_cor) {
            if (value_method == "max") {
                cor_merged[r, c] = max(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "min") {
                cor_merged[r, c] = min(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "mean") {
                cor_merged[r, c] = mean(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "range") {
                cor_merged[r, c] = max(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                }))) - min(unlist(lapply(cor_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else {
                stop("Invalid method for value merging")
            }
        }
    }
    return(cor_merged)
}

#' Merge clusters slot in domino object
#' 
#' Merge clusters based on selected method
#' 
#' @param dom_list A list of domino objects to be merged
#' @keywords internal
merge_clusters <- function(dom_list, feature_method) {
    clust_list = lapply(dom_list, slot, "clusters")
    clusters = lapply(clust_list, function(x) setNames(as.character(x), names(x)))
    clusters = factor(purrr::list_c(clusters))
    return(clusters)
}

#' Merge features slot of domino objects
#' 
#' Merge list of matrices of regulon scores by cell
#' 
#' @param dom_list A list of domino objects to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_features <- function(dom_list, feature_method) {
    feat_list = lapply(dom_list, slot, "features")
    # Return features based on feature_method with RVenn functions
    features = lapply(feat_list, rownames)
    features_list = RVenn::Venn(features)
    if (feature_method == "union") {
        features = RVenn::unite(features_list)
    } else if (feature_method == "intersect") {
        features = RVenn::overlap(features_list)
    } else if (feature_method == "outersect") {
        features = c()
        for (i in seq_along(feat_list)) {
            unique_items = RVenn::discern(features_list, i)
            features = c(features, unique_items)
        }
    } else {
        stop("Invalid method for feature merging")
    }
    subset = lapply(feat_list, function(x) {
        x[intersect(rownames(x), features), ]
    })
    # Then merge together by binding column-wise
    merged = do.call(cbind, subset)
    return(merged)
}

#' Merge complexes in linkages
#' 
#' Merge complexes, a list with names as complexes and values as genes
#' 
#' @param link_list A list of linkages to be merged
#' @keywords internal
merge_link_complex <- function(link_list) {
    # Get list of complexes from each item in link_list
    complex_list = lapply(link_list, function(x) {
        x["complexes"]$complexes
    })
    # Each item in list might have different complexes, so we need to get union of all complexes
    complexes = unique(unlist(lapply(complex_list, names)))
    # For each complex, get the genes in the complex
    complex_merged = list()
    for (c in complexes) {
        genes = unique(unlist(lapply(complex_list, function(x) {
            x[[c]]
        })))
        # Handle those semicolons, too:
        genes = vapply(genes, remove_semicolons, "character", USE.NAMES = FALSE)
        complex_merged[[c]] = genes
    }
    return(complex_merged)
}

#' Merge receptor ligand slot of linkages
#' 
#' Merge rec_lig slog, a list with receptor as names and ligands as values
#' 
#' @param link_list A list of linkages to be merged
#' @keywords internal
merge_link_rl <- function(link_list) {
    # Get list of rec_lig from each item in link_list
    rec_lig_list = lapply(link_list, function(x) {
        x["rec_lig"]$rec_lig
    })
    # Each item in list might have different receptors, so we need to get union of all receptors
    rec_lig = unique(unlist(lapply(rec_lig_list, names)))
    # For each receptor, get the ligands
    rec_lig_merged = list()
    for (r in rec_lig) {
        ligands = unique(unlist(lapply(rec_lig_list, function(x) {
            x[[r]]
        })))
        ligands = vapply(ligands, remove_semicolons, "character", USE.NAMES = FALSE)
        rec_lig_merged[[r]] = ligands
    }
    return(rec_lig_merged)
} 

#' Merge regulons in linkages
#' 
#' Merge tf_targets, a list with TF as names and targets as values
#' 
#' @param link_list A list of linkages to be merged
#' @keywords internal
merge_link_regulon <- function(link_list) {
    # Get list of tf_targets from each item in link_list
    tf_targets_list = lapply(link_list, function(x) {
        x["tf_targets"]$tf_targets
    })
    # Each item in list might have different TFs, so we need to get union of all TFs
    tf_targets = unique(unlist(lapply(tf_targets_list, names)))
    # For each TF, get the targets
    tf_targets_merged = list()
    for (tf in tf_targets) {
        ind = unlist(lapply(tf_targets_list, function(x) {
            which(tf %in% names(x))
        }))
        tf_targets_sub = tf_targets_list[names(ind)]
        targets = unique(unlist(lapply(tf_targets_sub, function(x) {
            x[[tf]]
        })))
        targets = vapply(targets, remove_semicolons, "character", USE.NAMES = FALSE)
        tf_targets_merged[[tf]] = targets
    }
    return(tf_targets_merged)
}

#' Merge tf_rec slot of linkages
#' 
#' Merge tf_rec, a list with TF as names and receptors as values
#' 
#' @param link_list A list of linkages to be merged
#' @keywords internal
merge_link_tfr <- function(link_list) {
    tf_rec_list = lapply(link_list, function(x) {
        x["tf_rec"]$tf_rec
    })
    tf_rec = unique(unlist(lapply(tf_rec_list, names)))
    tf_rec_merged = list()
    for (tf in tf_rec) {
        recs = unique(unlist(lapply(tf_rec_list, function(x) {
            x[[tf]]
        })))
        recs = vapply(recs, remove_semicolons, "character", USE.NAMES = FALSE)
        tf_rec_merged[[tf]] = recs
    }
    return(tf_rec_merged)
}

#' Merge clust_incoming_lig slot of linkages
#' 
#' Merge clust_incoming_lig, a list with cluster as names and ligands as values
#' 
#' @param link_list A list of linkages to be merged
#' @keywords internal
merge_link_in_lig <- function(link_list) {
    cl_lig_list = lapply(link_list, function(x) {
        x["clust_incoming_lig"]$clust_incoming_lig
    })
    cl_lig = unique(unlist(lapply(cl_lig_list, names)))
    cl_lig_merged = list()
    for (cl in cl_lig) {
        ligs = unique(unlist(lapply(cl_lig_list, function(x) {
            x[[cl]]
        })))
        ligs = vapply(ligs, remove_semicolons, "character", USE.NAMES = FALSE)
        cl_lig_merged[[cl]] = ligs
    }
    return(cl_lig_merged)
}

#' Merge clust_tf slot of linkages
#' 
#' Merge list with cluster as names and TFs as values
#' 
#' @param link_list A list of linkages to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_link_clust_tf <- function(link_list, feature_method) {
    cl_tf_list = lapply(link_list, function(x) {
        x["clust_tf"]$clust_tf
    })
    # Get cluster names (1 level down in list):
    clusts = unique(unlist(lapply(cl_tf_list, names)))
    # For each cluster, get elements from each list element
    cl_tf_merged = list()
    for (cl in clusts) {
        set = lapply(cl_tf_list, function(x) {
            x[[cl]]
        })
        set = RVenn::Venn(set)
        # Use RVenn to merge sets:
        if (feature_method == "union") {
            cl_tf_merged[[cl]] = RVenn::unite(set)
        } else if (feature_method == "intersect") {
            cl_tf_merged[[cl]] = RVenn::overlap(set)
        } else if (feature_method == "outersect") {
            cl_tf_merged[[cl]] = c()
            for (i in seq_along(cl_tf_list)) {
                unique_items = RVenn::discern(set, i)
                cl_tf_merged[[cl]] = c(cl_tf_merged[[cl]], unique_items)
            }
        } else {
            stop("Invalid method for feature merging")
        }
    }
    return(cl_tf_merged)
}

#' Merge clust_tf_rec slot of linkages
#' 
#' Merge list of clusters which is a list of TFs with receptors as elements
#' 
#' @param link_list A list of linkages to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_link_clust_tf_rec <- function(link_list, feature_method) {
    cl_tf_rec_list = lapply(link_list, function(x) {
        x["clust_tf_rec"]$clust_tf_rec
    })
    clusts = unique(unlist(lapply(cl_tf_rec_list, names)))
    cl_tf_rec_merged = list()
    for (cl in clusts) {
        # Grab the list of TFs and receptors for each cluster
        cl_tf_rec = lapply(cl_tf_rec_list, function(x) {
            names(x[[cl]])
        })
        # Then find which TFs to use based on feature_method using names
        tfs = RVenn::Venn(cl_tf_rec)
        if (feature_method == "intersect") {
            tfs_vec = RVenn::overlap(tfs)
        } else if (feature_method == "union") {
            tfs_vec = RVenn::unite(tfs)
        } else if (feature_method == "outersect") {
            tfs_vec = c()
            for (i in seq_along(cl_tf_rec)) {
                unique_items = RVenn::discern(tfs, i)
                tfs_vec = c(tfs_vec, unique_items)
            }
        } else {
            stop("Invalid method for feature merging")
        }
        # Then for each TF, select receptors based on feature_method
        cl_tf_rec_merged[[cl]] = list()
        if (!identical(tfs_vec, character(0))) {
            for (tf in tfs_vec) {
                recs = lapply(cl_tf_rec_list, function(x) {
                    x[[cl]][[tf]]
                })
                if (length(!is.null(recs)) > 1) {
                    recs = RVenn::Venn(recs)
                    if (feature_method == "intersect") {
                        recs_vec = RVenn::overlap(recs)
                    } else if (feature_method == "union") {
                        recs_vec = RVenn::unite(recs)
                    } else if (feature_method == "outersect") {
                        recs_vec = c()
                        for (i in seq_along(recs)) {
                            unique_items = RVenn::discern(recs, i)
                            recs_vec = c(recs_vec, unique_items)
                        }
                    } else {
                        stop("Invalid method for feature merging")
                    }
                } else {
                    recs_vec = unlist(recs[which(!is.null(recs))])
                }
                names(recs_vec) = NULL
                cl_tf_rec_merged[[cl]][[tf]] = recs_vec
            }
        } else {
            cl_tf_rec_merged[[cl]] = list()
        }
    }
    return(cl_tf_rec_merged)
}

#' Merge clust_rec slot of linkages
#' 
#' Merge list where name is cluster and items are receptors
#' 
#' @param link_list A list of linkages to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_link_clust_rec <- function(link_list, feature_method) {
    cl_rec_list = lapply(link_list, function(x) {
        x["clust_rec"]$clust_rec
    })
    clusts = unique(unlist(lapply(cl_rec_list, names)))
    cl_rec_merged = list()
    for (cl in clusts) {
        set = lapply(cl_rec_list, function(x) {
            x[[cl]]
        })
        set = RVenn::Venn(set)
        if (feature_method == "union") {
            cl_rec_merged[[cl]] = RVenn::unite(set)
        } else if (feature_method == "intersect") {
            cl_rec_merged[[cl]] = RVenn::overlap(set)
        } else if (feature_method == "outersect") {
            cl_rec_merged[[cl]] = c()
            for (i in seq_along(cl_rec_list)) {
                unique_items = RVenn::discern(set, i)
                cl_rec_merged[[cl]] = c(cl_rec_merged[[cl]], unique_items)
            }
        } else {
            stop("Invalid method for feature merging")
        }
    }
    return(cl_rec_merged)
}

#' Merge the linkages slot
#' 
#' Merge all linkages (complexes, rec_lig, tf_targets, clust_tf, tf_rec, clust_tf_rec, clust_rec, clust_incoming_lig)
#' 
#' @param dom_list A list of domino objects to be merged
#' @param feature_method A character string specifying the method for merging features
#' @keywords internal
merge_linkages <- function(dom_list, feature_method) {
    link_list = lapply(dom_list, slot, "linkages")
    merged_linkages = list(
        "complexes" = merge_link_complex(link_list),
        "rec_lig" = merge_link_rl(link_list),
        "tf_targets" = merge_link_regulon(link_list),
        "clust_tf" = merge_link_clust_tf(link_list, feature_method),
        "tf_rec" = merge_link_tfr(link_list),
        "clust_tf_rec" = merge_link_clust_tf_rec(link_list, feature_method),
        "clust_rec" = merge_link_clust_rec(link_list, feature_method),
        "clust_incoming_lig" = merge_link_in_lig(link_list)
        )
    return(merged_linkages)
}

#' Merge cluster signaling matrices
#' 
#' Merge lists with clusters as names, and matrix of receptors x ligand clusters as values
#' 
#' @param dom_list A list of domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_cl_signaling <- function(dom_list, value_method) {
    cl_signaling_list = lapply(dom_list, slot, "cl_signaling_matrices")
    clusts = unique(unlist(lapply(cl_signaling_list, names)))
    cl_signaling_merged = list()
    for (cl in clusts) {
        mat = lapply(cl_signaling_list, function(x) {
            x[[cl]]
        })
        # Build cl_signaling_merged[[cl]] empty matrix based no. of rows and columns
        row = unique(unlist(lapply(mat, rownames)))
        col = unique(unlist(lapply(mat, colnames)))
        cl_signaling_merged[[cl]] = matrix(0, nrow = length(row), ncol = length(col))
        rownames(cl_signaling_merged[[cl]]) = row
        colnames(cl_signaling_merged[[cl]]) = col
        # Build matrix for cluster based on value_method
        for (r in rownames(cl_signaling_merged[[cl]])) {
            for (c in colnames(cl_signaling_merged[[cl]])) {
                if (value_method == "max") {
                    cl_signaling_merged[[cl]][r, c] = max(unlist(lapply(mat, function(x) {
                        rc_val_mat(x, r, c)
                    })), na.rm = TRUE)
                } else if (value_method == "min") {
                    cl_signaling_merged[[cl]][r, c] = min(unlist(lapply(mat, function(x) {
                        rc_val_mat(x, r, c)
                    })), na.rm = TRUE)
                } else if (value_method == "mean") {
                    cl_signaling_merged[[cl]][r, c] = mean(unlist(lapply(mat, function(x) {
                        rc_val_mat(x, r, c)
                    })), na.rm = TRUE)
                } else if (value_method == "range") {
                    cl_signaling_merged[[cl]][r, c] = max(unlist(lapply(mat, function(x) {
                        rc_val_mat(x, r, c)
                    })), na.rm = TRUE) - min(unlist(lapply(sub, function(x) {
                        rc_val_mat(x, r, c)
                    })), na.rm = TRUE)
                } else {
                    stop("Invalid method for value merging")
                }
            }
        }
    }
    return(cl_signaling_merged)
}

#' Merge signaling matrix
#' 
#' Merge matrix of receiving clusters by sending clusters
#' 
#' @param dom_list A list of domino objects to be merged
#' @param value_method A character string specifying the method for merging numeric values
#' @keywords internal
merge_signaling <- function(dom_list, value_method) {
    signaling_list = lapply(dom_list, slot, "signaling")
    row = unique(unlist(lapply(signaling_list, rownames)))
    col = unique(unlist(lapply(signaling_list, colnames)))
    signaling_merged = matrix(0, nrow = length(row), ncol = length(col))
    rownames(signaling_merged) = row
    colnames(signaling_merged) = col
    for (r in rownames(signaling_merged)) {
        for (c in colnames(signaling_merged)) {
            if (value_method == "max") {
                signaling_merged[r, c] = max(unlist(lapply(signaling_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "min") {
                signaling_merged[r, c] = min(unlist(lapply(signaling_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "mean") {
                signaling_merged[r, c] = mean(unlist(lapply(signaling_list, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else if (value_method == "range") {
                signaling_merged[r, c] = max(unlist(lapply(signaling_list, function(x) {
                    rc_val_mat(x, r, c)
                }))) - min(unlist(lapply(sub, function(x) {
                    rc_val_mat(x, r, c)
                })))
            } else {
                stop("Invalid method for value merging")
            }
        }
    }
    return(signaling_merged)
}

#' Combine domino objects
#' Merge list of domino objects, with specified methods for feature and value merging
#'
#' @param dom_list A list of domino objects to be merged
#' @param feature_method A character string specifying the method for merging features
#'   (intersect, union, outersect)
#' @param value_method A character string specifying the method for merging numeric values
#'   (average, min, max, range)
#' @return A domino object with merged data across slots
#' @export
#' 
#' @examples 
#' dom_list = list(dom1, dom2, dom3)
#' combine_dom = merge_dom(dom_list, feature_method = "intersect", value_method = "mean")
#' contrast_dom = merge_dom(dom_list, feature_method = "outersect", value_method = "range")
merge_dom <- function(dom_list, feature_method = c("intersect", "union", "outersect"),
                        value_method = c("average", "min", "max", "range")) {
        dom = new("domino")
        slot(dom, "db_info") = merge_db_info(dom_list)
        slot(dom, "z_scores") = merge_zs(dom_list, feature_method)
        slot(dom, "counts") = merge_counts(dom_list, feature_method)
        slot(dom, "clusters") = merge_clusters(dom_list)
        slot(dom, "features") = merge_features(dom_list, feature_method)
        slot(dom, "cor") = merge_cor(dom_list, value_method)
        slot(dom, "linkages") = merge_linkages(dom_list, feature_method)
        slot(dom, "clust_de") = merge_de(dom_list, value_method)
        slot(dom, "misc") = merge_misc(dom_list, value_method)
        slot(dom, "cl_signaling_matrices") = merge_cl_signaling(dom_list, value_method)
        slot(dom, "signaling") = merge_signaling(dom_list, value_method)
        return(dom)
    }

#' Get value from matrix based on row and column index if present
#' 
#' Check whether row and column index exists in matrix; if so, pull value, and if
#' not, return 0 rather than throwing an error
#' 
#' @param matrix A matrix to pull values from
#' @param r_index A row index to pull value from
#' @param c_index A column index to pull value from
#' @return Value from matrix if present, 0 otherwise
#' @keywords internal
rc_val_mat <- function(matrix, r_index, c_index) {
    if (r_index %in% rownames(matrix) & c_index %in% colnames(matrix)) {
        return(matrix[r_index, c_index])
    } else {
        return(0)
    }
}