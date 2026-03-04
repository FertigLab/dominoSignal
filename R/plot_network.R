#' @import grid
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom igraph graph V E layout_in_circle layout_on_sphere layout_randomly layout_with_fr layout_with_kk simplify
#' @importFrom ggpubr ggscatter
#' @import grDevices
#'
NULL

#' Create a cluster to cluster signaling network diagram
#'
#' Creates a network diagram of signaling between clusters. Nodes are clusters
#' and directed edges indicate signaling from one cluster to another. Edges are
#' colored based on the color scheme of the ligand expressing cluster
#'
#' @param dom a domino object with network built ([build_domino()])
#' @param cols named vector indicating the colors for clusters. Values are colors and names must match clusters in
#'   the domino object. If left as NULL then ggplot colors are generated for the clusters
#' @param edge_weight weight for determining thickness of edges on plot. Signaling values are multiplied by this value
#' @param clusts vector of clusters to be included in the network plot
#' @param showOutgoingSignalingClusts vector of clusters to plot the outgoing signaling from
#' @param showIncomingSignalingClusts vector of clusters to plot the incoming signaling on
#' @param min_thresh minimum signaling threshold. Values lower than the threshold will be set to the threshold.
#'   Defaults to -Inf for no threshold
#' @param max_thresh maximum signaling threshold for plotting. Values higher than the threshold will be set to the
#'   threshold. Defaults to Inf for no threshold
#' @param normalize options to normalize the signaling matrix. Accepted inputs are 'none' for no normalization,
#'   'rec_norm' to normalize to the maximum value with each receptor cluster, or 'lig_norm' to normalize to the
#'   maximum value within each ligand cluster
#' @param scale how to scale the values (after thresholding). Options are 'none', 'sqrt' for square root,
#'   'log' for log10, or 'sq' for square
#' @param layout type of layout to use. Options are 'random', 'sphere', 'circle', 'fr' for Fruchterman-Reingold
#'   force directed layout, and 'kk' for Kamada Kawai for directed layout
#' @param scale_by how to size vertices. Options are 'lig_sig' for summed outgoing signaling, 'rec_sig' for summed
#'   incoming signaling, and 'none'. In the former two cases the values are scaled with asinh after summing all
#'   incoming or outgoing signaling
#' @param vert_scale integer used to scale size of vertices with our without variable scaling from size_verts_by.
#' @param plot_title text for the plot's title.
#' @param ... other parameters to be passed to plot when used with an igraph object.
#' @return An igraph plot rendered to the active graphics device
#' @export signaling_network
#' @examples
#' example(build_domino, echo = FALSE)
#' # basic usage
#' signaling_network(pbmc_dom_built_tiny, edge_weight = 2)
#' # scaling, thresholds, layouts, selecting clusters
#' signaling_network(
#'     pbmc_dom_built_tiny,
#'     showOutgoingSignalingClusts = "CD14_monocyte",
#'     scale = "none", normalize = "none", layout = "fr", scale_by = "none",
#'     vert_scale = 5, edge_weight = 2
#' )
#'
signaling_network <- function(
    dom, cols = NULL, edge_weight = 0.3, clusts = NULL, showOutgoingSignalingClusts = NULL,
    showIncomingSignalingClusts = NULL, min_thresh = -Inf, max_thresh = Inf, normalize = "none", scale = "sq",
    layout = "circle", scale_by = "rec_sig", vert_scale = 3, plot_title = NULL, ...
) {
    if (!length(dom@clusters)) {
        stop("This domino object was not built with clusters so there is no intercluster signaling.")
    }
    if (!dom@misc[["build"]]) {
        stop("Please build a signaling network with build_domino prior to plotting.")
    }
    # Get signaling matrix
    mat <- dom@signaling

    if (anyNA(mat)) {
        warning("Some values are NA, replacing with 0s.")
        mat[is.na(mat)] <- 0
    }

    if (!is.null(clusts)) {
        mat <- mat[paste0("R_", clusts), paste0("L_", clusts), drop = FALSE]
    }
    if (!is.null(showOutgoingSignalingClusts)) {
        mat <- mat[, paste0("L_", showOutgoingSignalingClusts), drop = FALSE]
    }
    if (!is.null(showIncomingSignalingClusts)) {
        mat <- mat[paste0("R_", showIncomingSignalingClusts), , drop = FALSE]
    }
    if (!any(mat > 0)) {
        warning("No signaling found")
        return(NULL)
    }
    if (is.null(cols)) {
        cols <- ggplot_col_gen(nlevels(dom@clusters))
        names(cols) <- levels(dom@clusters)
    }

    mat[which(mat > max_thresh)] <- max_thresh
    mat[which(mat < min_thresh)] <- min_thresh
    mat <- switch(scale,
        sqrt = sqrt(mat),
        log = log10(mat + 1),
        sq = mat^2,
        none = mat,
        stop("Do not recognize scale input")
    )
    mat <- switch(normalize,
        rec_norm = {
            if (ncol(mat) > 1) {
                do_norm(mat, "row")
            } else {
                mat
            }
        },
        lig_norm = {
            if (nrow(mat) > 1) {
                do_norm(mat, "col")
            } else {
                mat
            }
        },
        none = mat,
        stop("Do not recognize normalize input")
    )
    links <- character(0)
    weight <- numeric(0)
    for (rcl in rownames(mat)) {
        for (lcl in colnames(mat)) {
            if (mat[rcl, lcl] == 0) {
                next
            }
            lig_cl <- gsub("L_", "", lcl, fixed = TRUE)
            rec_cl <- gsub("R_", "", rcl, fixed = TRUE)
            links <- c(links, as.character(lig_cl), as.character(rec_cl))
            weight[paste0(lig_cl, "|", rec_cl)] <- mat[rcl, lcl]
        }
    }
    graph <- igraph::graph(links)
    # Get vert colors and scale size if desired.
    igraph::V(graph)$label.dist <- 1.5
    igraph::V(graph)$label.color <- "black"
    v_cols <- cols[names(igraph::V(graph))]
    if (scale_by == "lig_sig" && all(gsub("L_", "", colnames(mat), fixed = TRUE) %in% names(igraph::V(graph)))) {
        vals <- asinh(colSums(mat))
        vals <- vals[paste0("L_", names(igraph::V(graph)))]
        igraph::V(graph)$size <- vals * vert_scale
    } else if (scale_by == "rec_sig" && all(gsub("R_", "", rownames(mat), fixed = TRUE) %in% names(igraph::V(graph)))) {
        vals <- asinh(rowSums(mat))
        vals <- vals[paste0("R_", names(igraph::V(graph)))]
        igraph::V(graph)$size <- vals * vert_scale
    } else {
        igraph::V(graph)$size <- vert_scale
    }
    # Get vert angle for labeling circos plot
    if (layout == "circle") {
        v_angles <- seq_along(igraph::V(graph))
        v_angles <- -2 * pi * (v_angles - 1) / length(v_angles)
        igraph::V(graph)$label.degree <- v_angles
    }
    names(v_cols) <- NULL
    igraph::V(graph)$color <- v_cols
    # Get edge color. weights, and lines
    weight <- weight[attr(igraph::E(graph), "vnames")]
    e_cols <- character(0)
    for (e in names(weight)) {
        lcl <- strsplit(e, "|", fixed = TRUE)[[1]][1]
        e_cols <- c(e_cols, cols[lcl])
    }
    names(weight) <- NULL
    names(e_cols) <- NULL
    igraph::E(graph)$width <- weight * edge_weight
    igraph::E(graph)$color <- e_cols
    igraph::E(graph)$arrow.size <- 0
    igraph::E(graph)$curved <- 0.5
    # Get edge colors
    l <- switch(layout,
        random = igraph::layout_randomly(graph),
        circle = igraph::layout_in_circle(graph),
        sphere = igraph::layout_on_sphere(graph),
        fr = igraph::layout_with_fr(graph),
        kk = igraph::layout_with_kk(graph),
        stop("Do not recognize layout input")
    )
    plot(graph, layout = l, main = plot_title, ...)
}

#' Create a gene association network
#'
#' Create a gene association network for genes from a given cluster. The
#' selected cluster acts as the receptor for the gene association network, so
#' only ligands, receptors, and features associated with the receptor cluster
#' will be included in the plot.
#'
#' @param dom Domino object with network built ([build_domino()])
#' @param clust Receptor cluster to create the gene association network for. A vector of clusters may be provided.
#' @param OutgoingSignalingClust Vector of clusters to plot the outgoing signaling from
#' @param class_cols Named vector of colors used to color classes of vertices. Values must be colors and names must
#'   be classes ('rec', 'lig', and 'feat' for receptors, ligands, and features).
#' @param cols Named vector of colors for individual genes.
#'   Genes not included in this vector will be colored according to class_cols.
#' @param lig_scale FALSE or a numeric value to scale the size of ligand vertices based on
#'   z-scored expression in the data set.
#' @param layout Type of layout to use. Options are 'grid', 'random', 'sphere', 'circle', 'fr' for
#'   Fruchterman-Reingold force directed layout, and 'kk' for Kamada Kawai for directed layout.
#' @param ... Other parameters to pass to plot() with an [igraph](https://r.igraph.org/) object.
#'   See [igraph](https://r.igraph.org/) manual for options.
#' @return An igraph plot rendered to the active graphics device
#' @export gene_network
#' @examples
#' # basic usage
#' example(build_domino, echo = FALSE)
#' gene_network(
#'     pbmc_dom_built_tiny,
#'     clust = "CD8_T_cell",
#'     OutgoingSignalingClust = "CD14_monocyte"
#' )
#'
gene_network <- function(
    dom, clust = NULL, OutgoingSignalingClust = NULL,
    class_cols = c(lig = "#FF685F", rec = "#47a7ff", feat = "#39C740"),
    cols = NULL, lig_scale = 1, layout = "grid", ...
) {
    if (!dom@misc[["build"]]) {
        warning("Please build a signaling network with build_domino prior to plotting.")
    }
    if (!length(dom@clusters)) {
        warning("This domino object wasn't built with clusters. The global signaling network will be shown.")
        lig_scale <- FALSE
    }
    # Get connections between TF and recs for clusters
    if (length(dom@clusters)) {
        all_sums <- numeric(0)
        tfs <- character(0)
        cl_with_signaling <- character(0)
        for (cl in as.character(clust)) {
            # Check if signaling exists for target cluster
            mat <- dom@cl_signaling_matrices[[cl]]
            if (dim(mat)[1] == 0) {
                message("No signaling found for ", cl, " under build parameters.")
                next
            }
            all_sums <- c(all_sums, rowSums(mat))
            tfs <- c(tfs, dom@linkages$clust_tf[[cl]])
            cl_with_signaling <- c(cl_with_signaling, cl)
        }
        all_sums <- all_sums[!duplicated(names(all_sums))]
        # If no signaling for target clusters then don't do anything
        if (length(tfs) == 0) {
            message("No signaling found for provided clusters")
            return()
        }
    } else {
        tfs <- dom@linkages$clust_tf[["clust"]]
    }
    links <- character(0)
    all_recs <- character(0)
    all_tfs <- character(0)
    for (cl in as.character(clust)) {
        for (tf in tfs) {
            recs <- dom@linkages$clust_tf_rec[[cl]][[tf]]
            all_recs <- c(all_recs, recs)
            if (length(recs)) {
                all_tfs <- c(all_tfs, tf)
            }
            for (rec in recs) {
                links <- c(links, rec, tf)
            }
        }
    }
    all_recs <- unique(all_recs)
    all_tfs <- unique(all_tfs)
    # Recs to ligs
    if (length(dom@clusters)) {
        allowed_ligs <- character(0)
        for (cl in cl_with_signaling) {
            if (!is.null(OutgoingSignalingClust)) {
                OutgoingSignalingClust <- paste0("L_", OutgoingSignalingClust)
                mat <- dom@cl_signaling_matrices[[cl]][, OutgoingSignalingClust]
                if (is.null(dim(mat))) {
                    allowed_ligs <- names(mat[mat > 0])
                    all_sums <- mat[mat > 0]
                } else {
                    # Remove ligands with 0s for all clusters
                    allowed_ligs <- rownames(mat[rowSums(mat) > 0, ]) 
                    all_sums <- rowSums(mat[rowSums(mat) > 0, ])
                }
            } else {
                allowed_ligs <- rownames(dom@cl_signaling_matrices[[cl]])
            }
        }
    } else {
        allowed_ligs <- rownames(dom@z_scores)
    }
    # Remove ligs not expressed in data set if desired
    all_ligs <- character(0)
    for (rec in all_recs) {
        ligs <- dom@linkages$rec_lig[[rec]]
        for (lig in ligs) {
            if (any(allowed_ligs == lig)) {
                links <- c(links, lig, rec)
                all_ligs <- c(all_ligs, lig)
            }
        }
    }
    all_ligs <- unique(all_ligs)
    # Make the graph
    graph <- igraph::graph(links)
    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
    v_cols <- rep("#BBBBBB", length(igraph::V(graph)))
    names(v_cols) <- names(igraph::V(graph))
    v_cols[all_tfs] <- class_cols["feat"]
    v_cols[all_recs] <- class_cols["rec"]
    v_cols[all_ligs] <- class_cols["lig"]
    if (!is.null(cols)) {
        v_cols[names(cols)] <- cols
    }
    names(v_cols) <- NULL
    igraph::V(graph)$color <- v_cols
    v_size <- rep(10, length(igraph::V(graph)))
    names(v_size) <- names(igraph::V(graph))
    if (lig_scale) {
        all_sums <- all_sums[names(all_sums) %in% names(v_size)]
        v_size[names(all_sums)] <- 0.5 * all_sums * lig_scale
    }
    names(v_size) <- NULL
    igraph::V(graph)$size <- v_size
    igraph::V(graph)$label.degree <- pi
    igraph::V(graph)$label.offset <- 2
    igraph::V(graph)$label.color <- "black"
    igraph::V(graph)$frame.color <- "black"
    igraph::E(graph)$width <- 0.5
    igraph::E(graph)$arrow.size <- 0
    igraph::E(graph)$color <- "black"
    if (layout == "grid") {
        l <- matrix(0, ncol = 2, nrow = length(igraph::V(graph)))
        rownames(l) <- names(igraph::V(graph))
        l[all_ligs, 1] <- -0.75
        l[all_recs, 1] <- 0
        l[all_tfs, 1] <- 0.75
        l[all_ligs, 2] <- (seq_along(all_ligs) / mean(seq_along(all_ligs)) - 1) * 2
        l[all_recs, 2] <- (seq_along(all_recs) / mean(seq_along(all_recs)) - 1) * 2
        l[all_tfs, 2] <- (seq_along(all_tfs) / mean(seq_along(all_tfs)) - 1) * 2
        rownames(l) <- NULL
    } else {
        l <- switch(layout,
            random = igraph::layout_randomly(graph),
            circle = igraph::layout_in_circle(graph),
            sphere = igraph::layout_on_sphere(graph),
            fr = igraph::layout_with_fr(graph),
            kk = igraph::layout_with_kk(graph),
            stop("Do not recognize layout input")
        )
    }
    plot(graph, layout = l, main = paste0("Signaling ", OutgoingSignalingClust, " to ", clust), ...)
    return(invisible(list(graph = graph, layout = l)))
}
