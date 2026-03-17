#' Summarize linkages from multiple domino objects
#'
#' Creates a [linkage_summary()] object storing the linkages learned in different domino objects as nested
#'   lists to facilitate comparisons of networks learned by domino across subject covariates.
#'
#' @param domino_results list of domino result with one domino object per subject. Names from the list must
#'   match subject_names
#' @param subject_meta data frame that includes the subject features by which the objects could be grouped.
#'   The first column must be subject names
#' @param subject_names vector of subject names in domino_results. If NULL, defaults to first column of subject_meta.
#' @return A linkage summary class object consisting of nested lists of the active transcription factors,
#'   active receptors, and incoming ligands for each cluster across multiple domino results
#' @export
#' @examples
#' data(PBMC)
#' data(SCENIC)
#' data(CellPhoneDB)
#' data(DominoObjects)
#'
#' # create alternative clustering by shuffling cluster assignments
#' clusters_tiny_alt <- setNames(
#'     PBMC$clusters_tiny[c(121:240, 1:120, 241:360)],
#'     names(PBMC$clusters_tiny)
#' )
#' clusters_tiny_alt <- as.factor(clusters_tiny_alt)
#'
#' # build an alternative domino object
#' pbmc_dom_tiny_alt <- create_domino(
#'     rl_map = CellPhoneDB$rl_map_tiny,
#'     features = SCENIC$auc_tiny,
#'     counts = PBMC$count_tiny,
#'     z_scores = PBMC$zscore_tiny,
#'     clusters = clusters_tiny_alt,
#'     tf_targets = SCENIC$regulon_list_tiny,
#'     use_clusters = TRUE,
#'     use_complexes = TRUE,
#'     remove_rec_dropout = FALSE
#' )
#'
#' pbmc_dom_built_tiny_alt <- build_domino(
#'     dom = pbmc_dom_tiny_alt,
#'     min_tf_pval = .05,
#'     max_tf_per_clust = Inf,
#'     max_rec_per_tf = Inf,
#'     rec_tf_cor_threshold = .1,
#'     min_rec_percentage = 0.01
#' )
#'
#' # create a list of domino objects
#' dom_ls <- list(
#'     dom1 = DominoObjects$built_dom_tiny,
#'     dom2 = pbmc_dom_built_tiny_alt
#' )
#'
#' # compare the linkages across the two domino objects
#' meta_df <- data.frame("ID" = c("dom1", "dom2"), "group" = c("A", "B"))
#' summarize_linkages(
#'     domino_results = dom_ls, subject_meta = meta_df,
#'     subject_names = meta_df$ID
#' )
summarize_linkages <- function(domino_results, subject_meta, subject_names = NULL) {
    
    check_arg(domino_results, allow_class = "list", need_names = TRUE)
    check_arg(subject_meta, allow_class = "data.frame")
    
    if (!is(domino_results, "list")) {
        stop("domino_results must be provided as a named list where names correspond to subject names")
    }
    if (is.null(subject_names)) {
        subject_names <- subject_meta[, 1]
    }
    if (!any(subject_names %in% names(domino_results))) {
        stop("No provided subject names match names from the domino results list")
    }
    if (sum(!subject_names %in% names(domino_results))) {
        extra_names <- subject_names[!subject_names %in% names(domino_results)]
        warning("Provided subject names included names not present in domino_results: ", toString(extra_names))
        subject_names <- subject_names[subject_names %in% names(domino_results)]
    }
    if (length(subject_names) < length(names(domino_results))) {
        warning("Linkage summary includes results only for provided subject names: ", toString(subject_names))
        subject_meta <- subject_meta[subject_meta[, 1] %in% subject_names, ]
    }
    subject_linkages <- list()
    for (id in subject_names) {
        dom <- domino_results[[id]]
        clusters <- levels(dom@clusters)
        c_features <- list()
        for (cluster in clusters) {
            # list of t.factors active in cell type
            tfs <- dom@linkages$clust_tf[[cluster]]
            # obtain all receptors linked to active t. factors
            rec <- dom@linkages$clust_rec[[cluster]]
            # limit to unique entries
            tfs <- unique(tfs)
            rec <- unique(rec)
            # obtain all incoming ligands that interact with the receptors limited to those present in
            # data set
            lig <- rownames(dom@cl_signaling_matrices[[cluster]])
            # linkages of t.factors and receptors per cluster
            tfs_rec <- character(0)
            for (t in tfs) {
                for (r in dom@linkages$clust_tf_rec[[cluster]][[t]]) {
                    tfs_rec <- c(tfs_rec, t, r)
                }
            }
            rec_lig <- character(0)
            for (r in rec) {
                for (l in dom@linkages$rec_lig[[r]]) {
                    if (!l %in% lig) {
                        next
                    }
                    rec_lig <- c(rec_lig, r, l)
                }
            }
            # stitch linked t.factors-receptors, receptors-ligands
            int_tfs_rec <- character(0)
            int_rec_lig <- character(0)
            for (i in 0:((length(tfs_rec) / 2) - 1)) {
                # count by twos and paste together with a <- denoting direction
                s <- i * 2
                interact <- paste(tfs_rec[1 + s], tfs_rec[2 + s], sep = " <- ")
                int_tfs_rec <- c(int_tfs_rec, interact)
            }
            for (i in 0:((length(rec_lig) / 2) - 1)) {
                # count by twos and paste together with a '<-' denoting direction
                s <- i * 2
                interact <- paste(rec_lig[1 + s], rec_lig[2 + s], sep = " <- ")
                int_rec_lig <- c(int_rec_lig, interact)
            }
            # save the features of this cluster
            c_features[[cluster]] <- list(
                tfs = tfs, rec = rec, incoming_lig = lig, tfs_rec = int_tfs_rec,
                rec_lig = int_rec_lig
            )
        }
        subject_linkages[[id]] <- c_features
    }
    return(linkage_summary(
        subject_names = factor(subject_names),
        subject_meta = subject_meta,
        subject_linkages = subject_linkages))
}
