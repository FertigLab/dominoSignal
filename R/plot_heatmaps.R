#' @import grid
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom igraph graph V E layout_in_circle layout_on_sphere layout_randomly layout_with_fr layout_with_kk simplify
#' @importFrom ggpubr ggscatter
#' @import grDevices
#' 
NULL

#' Create a network heatmap
#'
#' Creates a heatmap of the signaling network. Alternatively, the network
#' matrix can be accessed directly in the signaling slot of a domino object using 
#' the [dom_signaling()] function.
#'
#' @param dom domino object with network built ([build_domino()])
#' @param clusts vector of clusters to be included. If NULL then all clusters are used.
#' @param min_thresh minimum signaling threshold for plotting. Defaults to -Inf for no threshold.
#' @param max_thresh maximum signaling threshold for plotting. Defaults to Inf for no threshold.
#' @param scale how to scale the values (after thresholding). Options are 'none', 'sqrt' for square root, or 'log' for log10.
#' @param normalize options to normalize the matrix. Normalization is done after thresholding and scaling. Accepted inputs are 'none' for no normalization, 'rec_norm' to normalize to the maximum value with each receptor cluster, or 'lig_norm' to normalize to the maximum value within each ligand cluster 
#' @param ... other parameters to pass to  [ComplexHeatmap::Heatmap()]
#' @return A heatmap rendered to the active graphics device
#' @export signaling_heatmap
#' @examples
#' example(build_domino, echo = FALSE)
#' #basic usage
#' signaling_heatmap(pbmc_dom_built_tiny)
#' #scale
#' signaling_heatmap(pbmc_dom_built_tiny, scale = "sqrt")
#' #normalize
#' signaling_heatmap(pbmc_dom_built_tiny, normalize = "rec_norm")
#'
signaling_heatmap <- function(
    dom, clusts = NULL, min_thresh = -Inf, max_thresh = Inf, scale = "none",
    normalize = "none", ...) {
  if (!dom@misc[["build"]]) {
    stop("Please run build_domino prior to generate signaling network.")
  }
  if (!length(dom@clusters)) {
    stop("This domino object wasn't built with clusters so intercluster signaling cannot be generated.")
  }
  mat <- dom@signaling

  if (!is.null(clusts)) {
    mat <- mat[paste0("R_", clusts), paste0("L_", clusts)]
  }
  mat[which(mat > max_thresh)] <- max_thresh
  mat[which(mat < min_thresh)] <- min_thresh
  if (scale == "sqrt") {
    mat <- sqrt(mat)
  } else if (scale == "log") {
    mat <- log10(mat)
  } else if (scale != "none") {
    stop("Do not recognize scale input")
  }
  if (normalize == "rec_norm") {
    mat <- do_norm(mat, "row")
  } else if (normalize == "lig_norm") {
    mat <- do_norm(mat, "col")
  } else if (normalize != "none") {
    stop("Do not recognize normalize input")
  }
  if (any(is.na(mat))) { 
    warning("Some values are NA, replacing with 0s.")
    mat[is.na(mat)] <- 0
  }

  Heatmap(
    mat,
    name = "collective\nsignaling",
    ...
  )
}

#' Create a cluster incoming signaling heatmap
#'
#' Creates a heatmap of a cluster incoming signaling matrix. Each cluster has a
#' list of ligands capable of activating its enriched transcription factors. The
#' function creates a heatmap of cluster average expression for all of those
#' ligands. A list of all cluster incoming signaling matrices can be found in
#' the cl_signaling_matrices slot of a domino option as an alternative to this
#' plotting function.
#'
#' @param dom Domino object with network built ([build_domino()])
#' @param rec_clust Which cluster to select as the receptor. Must match naming of clusters in the domino object.
#' @param clusts Vector of clusters to be included. If NULL then all clusters are used.
#' @param min_thresh Minimum signaling threshold for plotting. Defaults to -Inf for no threshold.
#' @param max_thresh Maximum signaling threshold for plotting. Defaults to Inf for no threshold.
#' @param scale How to scale the values (after thresholding). Options are 'none', 'sqrt' for square root, or 'log' for log10.
#' @param normalize Options to normalize the matrix. Accepted inputs are 'none' for no normalization, 'rec_norm' to normalize to the maximum value with each receptor cluster, or 'lig_norm' to normalize to the maximum value within each ligand cluster 
#' @param title Either a string to use as the title or a boolean describing whether to include a title. In order to pass the 'main' parameter to  [ComplexHeatmap::Heatmap()]  you must set title to FALSE.
#' @param ... Other parameters to pass to  [ComplexHeatmap::Heatmap()]. Note that to use the 'column_title' parameter of  [ComplexHeatmap::Heatmap()]  you must set title = FALSE
#' @return a Heatmap rendered to the active graphics device
#' @export incoming_signaling_heatmap
#' @examples
#' example(build_domino, echo = FALSE)
#' #incoming signaling of the CD8  T cells
#' incoming_signaling_heatmap(pbmc_dom_built_tiny, "CD8_T_cell")
#'
incoming_signaling_heatmap <- function(
    dom, rec_clust, clusts = NULL, min_thresh = -Inf, max_thresh = Inf,
    scale = "none", normalize = "none", title = TRUE, ...) {
  if (!dom@misc[["build"]]) {
    stop("Please run domino_build prior to generate signaling network.")
  }
  if (!length(dom@clusters)) {
    stop("This domino object wasn't build with clusters so cluster specific expression is not possible.")
  }
  mat <- dom@cl_signaling_matrices[[rec_clust]]
  if (dim(mat)[1] == 0) {
    message("No signaling found for this cluster under build parameters.")
    return()
  }

  if (!is.null(clusts)) {
    mat <- mat[, paste0("L_", clusts), drop = FALSE]
  }
  mat[which(mat > max_thresh)] <- max_thresh
  mat[which(mat < min_thresh)] <- min_thresh
  if (scale == "sqrt") {
    mat <- sqrt(mat)
  } else if (scale == "log") {
    mat <- log10(mat)
  } else if (scale != "none") {
    stop("Do not recognize scale input")
  }
  if (normalize == "rec_norm") {
    if (ncol(mat) > 1) {
      mat <- do_norm(mat, "row")
    }
  } else if (normalize == "lig_norm") {
    if (nrow(mat) > 1) {
      mat <- do_norm(mat, "col")
    }
  } else if (normalize != "none") {
    stop("Do not recognize normalize input")
  }

  if (any(is.na(mat))) { 
    warning("Some values are NA, replacing with 0s.")
    mat[is.na(mat)] <- 0
  }

  if (title == TRUE) {
    return(
      Heatmap(
        mat, 
        name = "expression",
        column_title = paste0("Expression of ligands targeting cluster ", rec_clust), 
        ...
      )
    )
  } else if (title == FALSE) {
    return(
      Heatmap(
        mat, 
        name = "expression",
        ...
      )
    )
  } else {
    return(
      Heatmap(
        mat,
        name = "expression",
        column_title = title, 
        ...
      )
    )
  }
}

#' Create a heatmap of features organized by cluster
#'
#' Creates a heatmap of transcription factor activation scores by cells grouped by cluster.
#'
#' @param dom Domino object with network built ([build_domino()])
#' @param bool Boolean indicating whether the heatmap should be continuous or boolean. If boolean then bool_thresh will be used to determine how to define activity as positive or negative.
#' @param bool_thresh Numeric indicating the threshold separating 'on' or 'off' for feature activity if making a boolean heatmap.
#' @param title Either a string to use as the title or a boolean describing whether to include a title. In order to pass the 'main' parameter to  [ComplexHeatmap::Heatmap()]  you must set title to FALSE.
#' @param norm Boolean indicating whether or not to normalize the transcrption factors to their max value.
#' @param feats Either a vector of features to include in the heatmap or 'all' for all features. If left NULL then the features selected for the signaling network will be shown.
#' @param ann_cols Boolean indicating whether to include cell cluster as a column annotation. Colors can be defined with cols. If FALSE then custom annotations can be passed to [ComplexHeatmap::Heatmap()].
#' @param cols Named vector of colors to annotate cells by cluster color. Values are taken as colors and names as cluster. If left as NULL then default ggplot colors will be generated.
#' @param min_thresh Minimum threshold for color scaling if not a boolean heatmap
#' @param max_thresh Maximum threshold for color scaling if not a boolean heatmap
#' @param ... Other parameters to pass to  [ComplexHeatmap::Heatmap()] . Note that to use the 'main' parameter of  [ComplexHeatmap::Heatmap()]  you must set title = FALSE and to use 'annCol' or 'annColors' ann_cols must be FALSE.
#' @return A heatmap rendered to the active graphics device
#' @export feat_heatmap
#' @examples 
#' #basic usage
#' example(build_domino, echo = FALSE)
#' feat_heatmap(pbmc_dom_built_tiny)
#' #using thresholds
#' feat_heatmap(
#'  pbmc_dom_built_tiny, min_thresh = 0.1, 
#'  max_thresh = 0.6, norm = TRUE, bool = FALSE)
#' 
feat_heatmap <- function(
    dom, feats = NULL, bool = FALSE, bool_thresh = 0.2, title = TRUE, norm = FALSE,
    cols = NULL, ann_cols = TRUE, min_thresh = NULL, max_thresh = NULL, ...) {
  if (!length(dom@clusters)) {
    warning("This domino object wasn't built with clusters. Cells will not be ordered.")
    ann_cols <- FALSE
  }
  mat <- dom@features
  cl <- dom@clusters
  cl <- sort(cl)
  if (norm & (!is.null(min_thresh) | !is.null(max_thresh))) {
    warning("You are using norm with min_thresh and max_thresh. Note that values will be thresholded AFTER normalization.")
  }
  if (norm) {
    mat <- do_norm(mat, "row")
  }
  if (!is.null(min_thresh)) {
    mat[which(mat < min_thresh)] <- min_thresh
  }
  if (!is.null(max_thresh)) {
    mat[which(mat > max_thresh)] <- max_thresh
  }
  if (bool) {
    cp <- mat
    cp[which(mat >= bool_thresh)] <- 1
    cp[which(mat < bool_thresh)] <- 0
    mat <- cp
  }
  if (title == TRUE) {
    title <- "Feature expression by cluster"
  }
  if (is.null(feats)) {
    feats <- c()
    links <- dom@linkages$clust_tf
    for (i in links) {
      feats <- c(feats, i)
    }
    feats <- unique(feats)
  } else if (feats[1] != "all") {
    mid <- match(feats, rownames(dom@features))
    na <- which(is.na(mid))
    na_feats <- paste(feats[na], collapse = " ")
    if (length(na) != 0) {
      message("Unable to find ", na_feats)
      feats <- feats[-na]
    }
  } else if (feats == "all") {
    feats <- rownames(mat)
  }
  if (length(cl)) {
    mat <- mat[feats, names(cl)]
  }
  if (ann_cols) {
    ac <- list(Cluster = cl)
    names(ac[[1]]) <- c()
    if (is.null(cols)) {
      cols <- ggplot_col_gen(length(levels(cl)))
      names(cols) <- levels(cl)
    }
    # cols <- list(Cluster = cols)
    feat_anno <- columnAnnotation(
      Cluster = cl,
      col = list(Cluster = cols)
    )
  }
  if (title != FALSE & ann_cols != FALSE) {
    Heatmap(
      mat,
      name = "feature\nactivity",
      top_annotation = feat_anno,
      cluster_columns = FALSE, show_column_names = FALSE,
      column_title = title,
      ...
    )
  } else if (title == FALSE & ann_cols != FALSE) {
    Heatmap(
      mat,
      name = "feature\nactivity",
      top_annotation = feat_anno,
      cluster_columns = FALSE, show_column_names = FALSE,
      ...
    )
  } else if (title != FALSE & ann_cols == FALSE) {
    Heatmap(
      mat,
      name = "feature\nactivity",
      top_annotation = feat_anno,
      cluster_columns = FALSE, show_column_names = FALSE,
      column_title = title,
      ...
    )
  } else if (title == FALSE & ann_cols == FALSE) {
    Heatmap(
      mat,
      name = "feature\nactivity",
      cluster_columns = FALSE, show_column_names = FALSE,
      ...
    )
  }
}

#' Create a heatmap of correlation between receptors and transcription factors
#'
#' Creates a heatmap of correlation values between receptors and transcription
#' factors either with boolean threshold or with continuous values displayed
#'
#' @param dom Domino object with network built ([build_domino()])
#' @param bool Boolean indicating whether the heatmap should be continuous or boolean. If boolean then bool_thresh will be used to determine how to define activity as positive or negative.
#' @param bool_thresh Numeric indicating the threshold separating 'on' or 'off' for feature activity if making a boolean heatmap.
#' @param title Either a string to use as the title or a boolean describing whether to include a title. In order to pass the 'main' parameter to  [ComplexHeatmap::Heatmap()]  you must set title to FALSE.
#' @param feats Either a vector of features to include in the heatmap or 'all' for all features. If left NULL then the features selected for the signaling network will be shown.
#' @param recs Either a vector of receptors to include in the heatmap or 'all' for all receptors. If left NULL then the receptors selected in the signaling network connected to the features plotted will be shown.
#' @param mark_connections Boolean indicating whether to add an 'x' in cells where there is a connected receptor or TF. Default FALSE.
#' @param ... Other parameters to pass to  [ComplexHeatmap::Heatmap()] . Note that to use the 'main' parameter of  [ComplexHeatmap::Heatmap()]  you must set title = FALSE and to use 'annCol' or 'annColors' ann_cols must be FALSE.
#' @return A heatmap rendered to the active graphics device
#' @export cor_heatmap
#' @examples 
#' example(build_domino, echo = FALSE)
#' #basic usage
#' cor_heatmap(pbmc_dom_built_tiny, title = "PBMC R-TF Correlations")
#' #show correlations above a specific value
#' cor_heatmap(pbmc_dom_built_tiny, bool = TRUE, bool_thresh = 0.1)
#' #identify combinations that are connected
#' cor_heatmap(pbmc_dom_built_tiny, bool = FALSE, mark_connections = TRUE)
#'  
cor_heatmap <- function(
    dom, bool = FALSE, bool_thresh = 0.15, title = TRUE, feats = NULL, recs = NULL,
    mark_connections = FALSE, ...) {
  mat <- dom@cor
  if (bool) {
    cp <- mat
    cp[which(mat >= bool_thresh)] <- 1
    cp[which(mat < bool_thresh)] <- 0
    mat <- cp
  }
  if (title == TRUE) {
    title <- "Correlation of features and receptors"
  }
  if (is.null(feats)) {
    feats <- c()
    links <- dom@linkages$clust_tf
    for (i in links) {
      feats <- c(feats, i)
    }
    feats <- unique(feats)
  } else if (feats[1] != "all") {
    mid <- match(feats, rownames(dom@features))
    na <- which(is.na(mid))
    na_feats <- paste(feats[na], collapse = " ")
    if (length(na) != 0) {
      message("Unable to find ", na_feats)
      feats <- feats[-na]
    }
  } else if (identical(feats, "all")) {
    feats <- rownames(mat)
  }
  if (is.null(recs)) {
    recs <- c()
    links <- dom@linkages$tf_rec
    for (feat in feats) {
      feat_recs <- links[[feat]]
      if (length(feat_recs) > 0) {
        recs <- c(recs, feat_recs)
      }
    }
    recs <- unique(recs)
  } else if (identical(recs, "all")) {
    recs <- rownames(mat)
  }
  mat <- mat[recs, feats]
  if (mark_connections) {
    cons <- mat
    cons[] <- ""
    for (feat in feats) {
      feat_recs <- dom@linkages$tf_rec[[feat]]
      if (length(feat_recs)) {
        cons[feat_recs, feat] <- "X"
      }
    }
  }
  if (title != FALSE & mark_connections) {
    Heatmap(
      mat,
      name = "rho",
      column_title = title,
      cell_fun = function(j, i, x, y, w, h, col){
        grid.text(
          cons[i,j], x, y,
          gp = gpar(col = "#000000")
        )
      },
      ...
    )
  } else {
    Heatmap(
      mat, 
      name = "rho",
      ...
    )
  }
}
