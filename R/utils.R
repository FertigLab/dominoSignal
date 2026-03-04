#' Check input arguments
#'
#' Accepts an object and rules to check against; stops if requirements are not met
#'
#' @param arg the argument to check
#' @param allow_class vector of allowed classes
#' @param allow_len vector of allowed lengths
#' @param allow_range range of minimum and maximum values i.e. c(1, 5)
#' @param allow_values vector of allowed values
#' @param need_vars vector of required variables
#' @param need_colnames vogical for whether colnames are required
#' @param need_rownames logical for whether rownames are required
#' @param need_names logical for whether names are required
#' @return Logical indicating whether the argument meets the requirements
#' @keywords internal
#'
check_arg <- function(arg, allow_class = NULL, allow_len = NULL,
    allow_range = NULL, allow_values = NULL, need_vars = NULL, need_colnames = FALSE,
    need_rownames = FALSE, need_names = FALSE) {

    argname <- deparse(substitute(arg))
    classes <- paste(allow_class, collapse = ",")
    lens <- paste(allow_len, collapse = ",")

    if (!is.null(allow_class) && !inherits(arg, allow_class, which = FALSE)) {
        stop(sprintf("Class of %s must be one of: %s", argname, classes))
    }

    if (!is.null(allow_len) && !(length(arg) %in% allow_len)) {
        stop(sprintf("Length of %s must be one of: %s", argname, lens))
    }

    if (!is.null(need_vars) && !all(need_vars %in% names(arg))) {
        stop(sprintf(
            "Required variables %s not found in %s",
            toString(need_vars), argname
        ))
    }

    if (need_rownames && is.null(rownames(arg))) {
        stop(sprintf("No rownames found in %s", argname))
    }

    if (need_colnames && is.null(colnames(arg))) {
        stop(sprintf("No colnames found in %s", argname))
    }

    if (need_names && is.null(names(arg))) {
        stop(sprintf("No names found in %s", argname))
    }

    if (!is.null(allow_range) && (any(arg < allow_range[1]) || any(arg > allow_range[2]))) {
        stop(sprintf(
            "All values in %s must be between %s and %s",
            argname, allow_range[1], allow_range[2]
        ))
    }

    if (!is.null(allow_values) && !all(arg %in% allow_values)) {
        stop(sprintf(
            "All values in %s must be one of: %s",
            argname, toString(allow_values)
        ))
    }
}

#' Read in data if an object looks like path to it
#'
#' @param obj object to read if not already object
#' @return Object itself or data read in from path
#' @keywords internal
read_if_char <- function(obj) {
    if (is(obj, "character")) {
        check_arg(obj, allow_class = "character", allow_len = 1)
        obj <- read.csv(obj, stringsAsFactors = FALSE)
    }
    return(obj)
}

#' Change cases of True/False syntax from Python to TRUE/FALSE R syntax
#'
#' @param obj object that will be converted
#' @return The converted object
#' @keywords internal
conv_py_bools <- function(obj) {
    for (x in colnames(obj)) {
        bools <- sort(unique(obj[[x]]))
        if (identical(bools, c("False", "True"))) {
            obj[[x]] <- obj[[x]] == "True"
        }
    }
    return(obj)
}

#' Convert between ligand names and gene names
#'
#' @param dom A domino object
#' @param genes A vector of genes on which to resolve ligand and gene names
#'
#' @return A vector of names where ligand names have been replaced with gene names if applicable
#' @keywords internal
resolve_names <- function(dom, genes) {
    rl_map <- dom@misc[["rl_map"]]
    genes_resolved <- vapply(genes, FUN.VALUE = character(1), FUN = function(l) {
        int <- rl_map[rl_map$L.name == l, ][1, ]
        if ((int$L.name != int$L.gene) & !grepl(",", int$L.gene, fixed = TRUE)) {
            int$L.gene
        } else {
            int$L.name
        }
    })
    return(genes_resolved)
}

#' Convert between complex names and gene names
#'
#' @param dom A domino object
#' @param genes A vector of genes, some of which may be complexes
#'
#' @return A list where any complexes are mapped to a vector of
#'   component genes. The list names are set to the input gene names.
#' @keywords internal
resolve_complexes <- function(dom, genes) {
    genes_list <- lapply(genes, function(l) {
        if (l %in% names(dom@linkages$complexes)) {
            return(dom@linkages$complexes[[l]])
        } else {
            return(l)
        }
    })
    names(genes_list) <- genes
    return(genes_list)
}

#' Pool items from list into vector
#'
#' Helper function to convert from a nested series of lists to a single vector.
#'
#' @param list List to pull items from
#' @param list_names Names of items in list to pool
#' @return A vector containing all items in the list by list_names
#' @keywords internal
#'
lc <- function(list, list_names) {
    vec <- NULL
    for (name in list_names) {
        vec <- c(vec, list[[name]])
    }
    return(vec)
}

#' Generate ggplot colors
#'
#' Accepts a number of colors to generate and generates a ggplot color spectrum.
#'
#' @param n Number of colors to generate
#' @return A vector of colors according to ggplot color generation.
#' @keywords internal
#'
ggplot_col_gen <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    return(grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)])
}
