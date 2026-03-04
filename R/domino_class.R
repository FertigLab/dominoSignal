#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#'
NULL
#' The domino class
#'
#' The domino class contains all information necessary to calculate receptor-ligand
#' signaling. It contains z-scored expression, cell cluster labels, feature values,
#' and a referenced receptor-ligand database formatted as a receptor-ligand map.
#' Calculated intermediate values are also stored.
#'
#' @slot db_info List of data sets from ligand - receptor database
#' @slot counts Raw count gene expression data
#' @slot z_scores Matrix of z-scored expression data with cells as columns
#' @slot clusters Named factor with cluster identity of each cell
#' @slot features Matrix of features (TFs) to correlate receptor - ligand expression with. Cells are columns and features are rows.
#' @slot cor Correlation matrix of receptor expression to features.
#' @slot linkages List of lists containing info linking cluster->tf->rec->lig
#' @slot clust_de Data frame containing differential expression results for features by cluster.
#' @slot misc List of miscellaneous info pertaining to run parameters etc.
#' @slot cl_signaling_matrices Incoming signaling matrix for each cluster
#' @slot signaling Signaling matrix between all clusters.
#' @name domino-class
#' @rdname domino-class
#' @exportClass domino
#' @return An instance of class `domino `
#'
domino <- methods::setClass(
  Class = "domino",
  slots = c(
    db_info = "list",
    z_scores = "matrix",
    counts = "dgCMatrix",
    clusters = "factor",
    features = "matrix",
    cor = "matrix",
    linkages = "list",
    clust_de = "matrix",
    misc = "list",
    cl_signaling_matrices = "list",
    signaling = "matrix"
  ),
  prototype = list(
    misc = list("build" = FALSE)
  )
)

#' Print domino object
#'
#' Prints a summary of a domino object
#'
#' @param x A domino object
#' @param ... Additional arguments to be passed to other methods
#' @return A printed description of the number of cells and clusters in the domino object
#' @export
#' @examples
#' example(build_domino, echo = FALSE)
#' print(pbmc_dom_built_tiny)
#'
setMethod("print", "domino", function(x, ...) {
  if (x@misc$build) {
    message(
      "A domino object of ", length(x@clusters), " cells
                Contains signaling between ",
      length(levels(x@clusters)), " clusters
                Built with a maximum of ", x@misc$build_vars["max_tf_per_clust"],
      " TFs per cluster
                and a maximum of ", x@misc$build_vars["max_rec_per_tf"],
      " receptors per TF\n"
    )
  } else {
    message(c("A domino object of ", length(x@clusters), " cells\n", "A signaling network has not been built\n"),
      sep = ""
    )
  }
})
#' Show domino object information
#'
#' Shows content overview of domino object
#'
#' @param object A domino object
#' @return A printed description of cell numbers and clusters in the object
#' @export
#' @examples
#' example(build_domino, echo = FALSE)
#' show(pbmc_dom_built_tiny)
#' 
setMethod("show", "domino", function(object) {
  if (object@misc$build) {
    cat(c(
      "A domino object of ", length(object@clusters), " cells\n", "Built with signaling between ",
      length(levels(object@clusters)), " clusters\n"
    ), sep = "")
  } else {
    cat(c("A domino object of ", length(object@clusters), " cells\n", "A signaling network has not been built\n"),
      sep = ""
    )
  }
})
