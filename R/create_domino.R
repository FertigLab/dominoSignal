#' Create a domino object and prepare it for network construction
#'
#' This function reads in a receptor ligand signaling database, cell level
#' features of some kind (ie. output from pySCENIC), z-scored single cell data,
#' and cluster id for single cell data, calculates a correlation matrix between
#' receptors and other features (this is transcription factor module scores if
#' using pySCENIC), and finds features enriched by cluster. It will return a
#' domino object prepared for [build_domino()], which will calculate a signaling
#' network.
#'
#' @param rl_map Data frame where each row describes a receptor-ligand interaction with required columns
#'   gene_A & gene_B including the gene names for the receptor and ligand and type_A & type_B annotating
#'   if genes A and B are a ligand (L) or receptor (R)
#' @param features Either a path to a csv containing cell level features of interest (ie. the auc matrix from pySCENIC)
#'   or named matrix with cells as columns and features as rows.
#' @param counts Counts matrix for the data. This is only used to threshold receptors on dropout.
#' @param z_scores Matrix containing z-scored expression data for all cells with cells as columns and features as rows.
#' @param clusters Named factor containing cell cluster with names as cells.
#' @param use_clusters Boolean indicating whether to use clusters.
#' @param tf_targets Optional. A list where names are transcription factors and the stored values are character vectors
#'   of genes in the transcription factor's regulon.
#' @param verbose Boolean indicating whether or not to print progress during computation.
#' @param use_complexes Boolean indicating whether you wish to use receptor/ligand complexes in the receptor ligand
#'   signaling database. If FALSE, receptor/ligand pairs where either functions as a protein complex will not be
#'   considered when constructing the signaling network.
#' @param rec_min_thresh Minimum expression level of receptors by cell. Default is 0.025 or 2.5 percent of all cells
#'   in the data set. This is important when calculating correlation to connect receptors to transcription activation.
#'   If this threshold is too low then correlation calculations will proceed with very few cells with non-zero
#'   expression.
#' @param remove_rec_dropout Whether to remove receptors with 0 expression counts when calculating correlations.
#'   This can reduce false positive correlation calculations when receptors have high dropout rates.
#' @param tf_selection_method Selection of which method to target transcription factors. If 'clusters' then
#'   differential expression for clusters will be calculated. If 'variable' then the most variable transcription factors
#'   will be selected. If 'all' then all transcription factors in the feature matrix will be used. Default is
#'   'clusters'. Note that if you wish to use clusters for intercellular signaling downstream to MUST choose clusters.
#' @param tf_variance_quantile What proportion of variable features to take if using variance to threshold features.
#'   Default is 0.5. Higher numbers will keep more features. Ignored if tf_selection_method is not 'variable'
#' @return A domino object
#' @export create_domino
#' @examples
#' example(create_rl_map_cellphonedb, echo = FALSE)
#' example(create_regulon_list_scenic, echo = FALSE)
#' data(SCENIC)
#' data(PBMC)
#'
#' pbmc_dom_tiny <- create_domino(
#'  rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
#'  counts = PBMC$RNA_count_tiny, z_scores = PBMC$RNA_zscore_tiny,
#'  clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
#'  use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE,
#'  verbose = FALSE
#'  )
#'
#' pbmc_dom_tiny_no_clusters <- create_domino(
#'  rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
#'  counts = PBMC$RNA_count_tiny, z_scores =PBMC$RNA_zscore_tiny,
#'  clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
#'  use_clusters = FALSE, use_complexes = FALSE,
#'  rec_min_thresh = 0.1, remove_rec_dropout = TRUE,
#'  tf_selection_method = "all",
#'  verbose = FALSE
#'  )
#'
create_domino <- function(
    rl_map, features, counts = NULL, z_scores = NULL,
    clusters = NULL, use_clusters = TRUE, tf_targets = NULL, verbose = TRUE,
    use_complexes = TRUE, rec_min_thresh = 0.025, remove_rec_dropout = TRUE,
    tf_selection_method = "clusters", tf_variance_quantile = 0.5
) {
    # Check inputs:
    check_arg(rl_map,
        allow_class = "data.frame",
        need_vars = c("gene_A", "gene_B", "type_A", "type_B")
    )

    check_arg(features, allow_class = c("data.frame", "character", "matrix"))
    if (is.data.frame(features) || is.matrix(features)) {
        check_arg(features, need_rownames = TRUE, need_colnames = TRUE)
    }


    check_arg(counts,
        allow_class = c("matrix", "data.frame", "Matrix", "dgCMatrix"),
        need_rownames = TRUE, need_colnames = TRUE
    )
    check_arg(z_scores,
        allow_class = "matrix", need_rownames = TRUE,
        need_colnames = TRUE
    )
    check_arg(clusters, allow_class = "factor", need_names = TRUE)


    check_arg(rec_min_thresh, allow_class = "numeric", allow_range = c(0, 1))

    check_arg(tf_selection_method,
        allow_values = c("clusters", "variable", "all")
    )

    # Create object
    dom <- domino()
    dom@misc[["create"]] <- TRUE
    dom@misc[["build"]] <- FALSE
    dom@misc[["build_vars"]] <- NULL

    # Read in lr db info
    if (verbose) {
        message("Reading in and processing signaling database")
    }
    if ("database_name" %in% colnames(rl_map)) {
        dom@db_info <- rl_map
        if (verbose) {
            message("Database provided from source: ", unique(rl_map[["database_name"]]))
        }
    } else {
        dom@db_info <- rl_map
    }
    # check for receptors that match receptor complex syntax of comma seperated genes
    non_complex_index <- which(!grepl(",", rl_map[["gene_A"]], fixed = TRUE) &
            !grepl(",", rl_map[["gene_B"]], fixed = TRUE))
    # discard interactions including complexes if requested
    if (!use_complexes) {
        rl_map <- rl_map[non_complex_index, ]
    }
    # Get genes for receptors
    rl_reading <- NULL
    for (i in seq_len(nrow(rl_map))) {
        rl <- list()
        inter <- rl_map[i, ]
        ps <- ifelse(inter[["type_A"]] == "R", "A", "B")
        qs <- ifelse(ps == "A", "B", "A")
        R.gene <- inter[[paste0("gene_", ps)]]
        L.gene <- inter[[paste0("gene_", qs)]]
        rl[["R.gene"]] <- R.gene
        rl[["L.gene"]] <- L.gene
        if (paste0("uniprot_", ps) %in% names(inter)) {
            rl[["R.uniprot"]] <- inter[[paste0("uniprot_", ps)]]
        }
        if (paste0("uniprot_", qs) %in% names(inter)) {
            rl[["L.uniprot"]] <- inter[[paste0("uniprot_", qs)]]
        }
        if (paste0("name_", ps) %in% names(inter)) {
            rl[["R.name"]] <- inter[[paste0("name_", ps)]]
        }
        if (paste0("name_", qs) %in% names(inter)) {
            rl[["L.name"]] <- inter[[paste0("name_", qs)]]
        }
        rl <- as.data.frame(rl)
        rl_reading <- rbind(rl_reading, rl)
    }
    if (nrow(rl_reading) == 0) stop("No genes annotated as receptors included in rl_map")
    # save a list of complexes and their components
    dom@linkages$complexes <- NULL
    if (use_complexes) {
        complex_list <- list()
        for (i in seq_len(nrow(rl_reading))) {
            inter <- rl_reading[i, ]
            if (grepl(",", inter[["L.gene"]], fixed = TRUE)) {
                complex_list[[inter[["L.name"]]]] <- unlist(strsplit(inter[["L.gene"]], split = ",", fixed = TRUE))
            }
            if (grepl(",", inter[["R.gene"]], fixed = TRUE)) {
                complex_list[[inter[["R.name"]]]] <- unlist(strsplit(inter[["R.gene"]], split = ",", fixed = TRUE))
            }
        }
        dom@linkages$complexes <- complex_list
    }
    rec_genes <- unique(unlist(strsplit(rl_reading[["R.gene"]], split = ",", fixed = TRUE)))
    rec_names <- rl_reading[["R.name"]]
    # building RL linkages
    rec_lig_linkage <- list()
    for (rec in rec_names) {
        inter <- rl_reading[rl_reading[["R.name"]] == rec, ]
        ligs <- inter[["L.name"]]
        rec_lig_linkage[[rec]] <- ligs
    }
    dom@linkages[["rec_lig"]] <- rec_lig_linkage
    dom@misc[["rl_map"]] <- rl_reading
    # Get z-score and cluster info
    if (verbose) {
        message("Getting z_scores, clusters, and counts")
    }
    dom@z_scores <- z_scores
    if (!is.null(clusters)) {
        dom@clusters <- clusters
    }
    # Read in features matrix and calculate differential expression by cluster.
    if (is(features, "character")) {
        features <- read.csv(features, row.names = 1, check.names = FALSE)
    }
    features <- features[, colnames(dom@z_scores)]
    dom@features <- as.matrix(features)
    if (tf_selection_method == "clusters") {
        p_vals <- matrix(1, nrow = nrow(features), ncol = nlevels(dom@clusters))
        rownames(p_vals) <- rownames(features)
        colnames(p_vals) <- levels(dom@clusters)
        if (verbose) {
            message("Calculating feature enrichment by cluster")
            clust_n <- nlevels(dom@clusters)
        }
        for (clust in levels(dom@clusters)) {
            if (verbose) {
                cur <- which(levels(dom@clusters) == clust)
                message(cur, " of ", clust_n)
            }
            cells <- which(dom@clusters == clust)
            for (feat in rownames(dom@features)) {
                p_vals[feat, clust] <- stats::wilcox.test(
                    dom@features[feat, cells], dom@features[feat, -cells],
                    alternative = "g"
                )$p.value
            }
        }
        dom@clust_de <- p_vals
    }
    if (tf_selection_method == "all") {
        dom@clusters <- factor()
    }
    if (tf_selection_method == "variable") {
        dom@clusters <- factor()
        variances <- apply(dom@features, 1, function(x) {
            sd(x) / mean(x)
        })
        keep_n <- length(variances) * tf_variance_quantile
        keep_id <- which(rank(variances) > keep_n)
        dom@features <- dom@features[names(keep_id), ]
    }
    # store tf_targets in linkages if they are provided as a list
    if (is(tf_targets, "list")) {
        dom@linkages[["tf_targets"]] <- tf_targets
    } else {
        dom@linkages[["tf_targets"]] <- NULL
        message("tf_targets is not a list. No regulons stored")
    }
    # Calculate correlation matrix between features and receptors.
    dom@counts <- counts
    zero_sum <- rowSums(counts == 0)
    keeps <- which(zero_sum < (1 - rec_min_thresh) * ncol(counts))
    ser_receptors <- intersect(names(keeps), rec_genes)
    rho <- matrix(0, nrow = length(ser_receptors), ncol = nrow(dom@features))
    rownames(rho) <- ser_receptors
    colnames(rho) <- rownames(dom@features)
    if (verbose) {
        message("Calculating correlations")
        n_tf <- nrow(dom@features)
    }
    for (module in rownames(dom@features)) {
        # If df is provided then check if receptors are targets of TF. If they are then set
        # correlation equal to 0.
        if (verbose) {
            cur <- which(rownames(dom@features) == module)
            message(cur, " of ", n_tf)
        }
        if (!is.null(dom@linkages$tf_targets)) {
            tf <- gsub(pattern = "...", replacement = "", module, fixed = TRUE)
            # correction for AUC values from pySCENIC that append an elipses to TF names due to (+) characters
            # in the orignial python output
            module_targets <- tf_targets[[tf]]
            module_rec_targets <- intersect(module_targets, ser_receptors)
        } else {
            module_rec_targets <- NULL
        }
        scores <- dom@features[module, ]
        rhorow <- rep(0, length(ser_receptors))
        names(rhorow) <- ser_receptors
        for (rec in ser_receptors) {
            if (remove_rec_dropout) {
                keep_id <- which(dom@counts[rec, ] > 0)
                rec_z_scores <- dom@z_scores[rec, keep_id]
                tar_tf_scores <- scores[keep_id]
            } else {
                rec_z_scores <- dom@z_scores[rec, ]
                tar_tf_scores <- scores
            }
            # There are some cases where all the tfs are zero for the cells left after trimming
            # dropout for receptors. Skip those and set cor to zero manually.
            if (sum(tar_tf_scores) == 0) {
                rhorow[rec] <- 0
                next
            }
            corr <- stats::cor.test(
                rec_z_scores, tar_tf_scores,
                method = "spearman",
                alternative = "greater", exact = FALSE
            )
            rhorow[rec] <- corr$estimate
        }
        if (length(module_rec_targets) > 0) {
            rhorow[module_rec_targets] <- 0
        }
        rho[, module] <- rhorow
    }
    colnames(rho) <- rownames(dom@features)
    dom@misc$rec_cor <- rho
    # assess correlation among genes in the same receptor complex
    cor_list <- list()
    for (i in seq_along(names(dom@linkages$rec_lig))) {
        r <- names(dom@linkages$rec_lig)[i]
        if (r %in% names(dom@linkages$complexes)) {
            r_genes <- dom@linkages$complexes[[r]]
        } else {
            r_genes <- r
        }
        if (sum(rownames(rho) %in% r_genes) != length(r_genes)) {
            cor_list[[r]] <- rep(0, ncol(rho))
            next
        }
        if (length(r_genes) > 1) {
            gene_cor <- rho[rownames(rho) %in% r_genes, ]
            cor_med <- apply(gene_cor, 2, median)
            cor_list[[r]] <- cor_med
        } else {
            cor_list[[r]] <- rho[rownames(rho) == r_genes, ]
        }
    }
    c_cor <- t(as.data.frame(cor_list))
    dom@cor <- c_cor
    # If cluster methods are used, calculate percentage of non-zero expression of receptor genes
    # in clusters
    if (tf_selection_method == "clusters") {
        cl_rec_percent <- NULL
        for (rec in ser_receptors) {
            rec_percent <- vapply(X = levels(dom@clusters), FUN.VALUE = numeric(1), FUN = function(x) {
                # percentage of cells in cluster with non-zero expression of receptor gene
                sum(dom@counts[rec, dom@clusters == x] > 0) / length(dom@counts[rec, dom@clusters == x])
            })
            cl_rec_percent <- rbind(cl_rec_percent, rec_percent)
        }
        rownames(cl_rec_percent) <- ser_receptors
        dom@misc$cl_rec_percent <- cl_rec_percent
    }
    return(dom)
}
