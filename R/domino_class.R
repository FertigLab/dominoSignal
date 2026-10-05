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
#' @slot features Matrix of features (TFs) to correlate receptor - ligand expression with.
#'   Cells are columns and features are rows.
#' @slot cor Correlation matrix of receptor expression to features.
#' @slot linkages List of lists containing info linking cluster->tf->rec->lig
#' @slot clust_de Data frame containing differential expression results for features by cluster.
#' @slot misc List of miscellaneous info pertaining to run parameters etc., including the
#'   dominoSignal versions used to create (`create_version`) and build (`build_version`) the object.
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

valid_domino_misc <- function(misc) {
    # Run state flags must be single TRUE/FALSE values when present
    flags <- c("create", "build")
    bad_flag <- vapply(flags, function(f) {
        !is.null(misc[[f]]) && !isTRUE(misc[[f]]) && !isFALSE(misc[[f]])
    }, logical(1))
    err <- sprintf("misc$%s must be a single TRUE or FALSE", flags[bad_flag])
    # Package versions are optional (objects made before they were recorded lack them)
    versions <- c("create_version", "build_version")
    bad_version <- vapply(versions, function(v) {
        val <- misc[[v]]
        !is.null(val) && !(is.character(val) && length(val) == 1 && !is.na(package_version(val, strict = FALSE)))
    }, logical(1))
    c(err, sprintf("misc$%s must be a single valid version string", versions[bad_version]))
}

valid_domino_slots <- function(object) {
    # Cell-level matrices must have the same cells in the same order as clusters, since cells are indexed by position
    cells <- names(object@clusters)
    cell_slots <- c("z_scores", "counts", "features")
    bad_cells <- vapply(cell_slots, function(s) {
        m <- slot(object, s)
        length(cells) > 0 && ncol(m) > 0 && !identical(colnames(m), cells)
    }, logical(1))
    err <- sprintf("%s columns must be the same cells in the same order as names(clusters)", cell_slots[bad_cells])
    if (ncol(object@clust_de) > 0 && !setequal(colnames(object@clust_de), levels(object@clusters))) {
        err <- c(err, "clust_de columns must match levels(clusters)")
    }
    if (ncol(object@cor) > 0 && !setequal(colnames(object@cor), rownames(object@features))) {
        err <- c(err, "cor columns must match features rows")
    }
    err
}

valid_domino <- function(object) {
    err <- c(valid_domino_misc(object@misc), valid_domino_slots(object))
    if (length(err) == 0) TRUE else err
}

setValidity("domino", valid_domino)

#' Print domino object
#'
#' Prints a summary of a domino object
#'
#' @param x A domino object
#' @param ... Additional arguments to be passed to other methods
#' @return A printed description of the number of cells and clusters in the domino object
#' @export
#' @examples
#' data(DominoObjects)
#' print(DominoObjects$built_dom_tiny)
#'
setMethod("print", "domino", function(x, ...) {
    version_msg <- paste0(
        if (!is.null(x@misc$create_version)) paste0("Created with dominoSignal v", x@misc$create_version, "\n"),
        if (!is.null(x@misc$build_version)) paste0("Built with dominoSignal v", x@misc$build_version, "\n")
    )
    if (x@misc$build) {
        message(
            "A domino object of ", length(x@clusters), " cells
                Contains signaling between ",
            nlevels(x@clusters), " clusters
                Built with a maximum of ", x@misc$build_vars["max_tf_per_clust"],
            " TFs per cluster
                and a maximum of ", x@misc$build_vars["max_rec_per_tf"],
            " receptors per TF\n", version_msg
        )
    } else {
        message("A domino object of ", length(x@clusters), " cells\n", "A signaling network has not been built\n",
            version_msg
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
#' data(DominoObjects)
#' show(DominoObjects$built_dom_tiny)
#'
setMethod("show", "domino", function(object) {
    print(object)
})
