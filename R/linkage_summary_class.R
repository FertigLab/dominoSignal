#' The domino linkage summary class
#'
#' The linkage summary class contains linkages established in multiple domino
#' objects through gene regulatory network inference and reference to receptor-
#' ligand databases. A data frame summarizing meta features that describe the
#' domino objects compared in the linkage summary facilitates comparisons of
#' established linkages and differential signaling interactions across categorical
#' sample covariates.
#'
#' @slot subject_names unique names for each domino result included in the summary
#' @slot subject_meta data.frame with each row describing one subject and columns describing features of the
#'   subjects by which to draw comparisons of signaling networks
#' @slot subject_linkages nested list of linkages inferred for each subject. Lists are stored in a
#'   hierarchical structure of subject-cluster-linkage where linkages include transcription factors (tfs)
#'   linkages between transcription factors and receptors (tfs_rec), active receptors (rec), possible
#'   receptor-ligand interactions (rec_lig), and incoming ligands (incoming_lig)
#' @name linkage_summary-class
#' @rdname linkage_summary-class
#' @exportClass linkage_summary
#' @return an instance of class `linkage_summary`
#'
linkage_summary <- setClass(
    Class = "linkage_summary",
    slots = c(
        subject_names = "factor",
        subject_meta = "data.frame",
        subject_linkages = "list"
    )
)

#' Print linkage summary object
#'
#' Prints a description of a linkage summary object
#'
#' @param x A linkage summary object
#' @param ... Additional arguments to be passed to other methods
#' @return A printed description of the number of subjects, groups, and clusters in the linkage summary object
#' @export
#' @examples
#' data(LinkageSummary)
#' print(LinkageSummary$linkage_sum_tiny)
#'
setMethod("print", "linkage_summary", function(x, ...) {
    message("A linkage summary object of ", nlevels(slot(x, "subject_names")), " subjects with ",
        ncol(slot(x, "subject_meta")), " metadata annotations and linkages between ",
        max(lengths(slot(x, "subject_linkages"))), " clusters.")
})

#' Show linkage_summary object information
#'
#' Shows content overview of linkage_summary object
#'
#' @param object A linkage_summary object
#' @return A printed description of the number of subjects, groups, and clusters in the linkage summary object
#' @export
#' @examples
#' data(LinkageSummary)
#' LinkageSummary$linkage_sum_tiny
#'
setMethod("show", "linkage_summary", function(object) {
    message("A linkage summary object of ", nlevels(slot(object, "subject_names")), " subjects with ",
        ncol(slot(object, "subject_meta")), " metadata annotations and linkages between ",
        max(lengths(slot(object, "subject_linkages"))), " clusters.")
})

#' Subset a linkage_summary object
#' 
#' Subsets a linkage summary object by subject names or metadata
#' 
#' @param x A linkage_summary object
#' @param subset A logical expression referencing `subject_names` or columns of `subject_meta`
#' @return A linkage_summary object containing only the subjects that match the specified criteria.
#' @export
#' @examples 
#' data(LinkageSummary)
#' links <- LinkageSummary$linkage_sum_tiny
#' subset(links, subset = !subject_names %in% c("P1", "P2"))
#' subset(links, subset = group == "G1")
#' subset(links, subset = subject_names %in% c("P1", "P2", "P4") | group == "G2")

setMethod("subset", "linkage_summary", function(x, subset) {
    check_arg(x, allow_class = "linkage_summary", allow_len = 1)
    # Capture the subset argument unevaluated; it references subject_names/subject_meta
    # columns rather than variables in the calling environment
    subset_call <- substitute(subset)
    check_arg(subset_call, allow_class = "call")
    # Only allow the subset expression to reference subject_names or columns of subject_meta
    valid_vars <- c("subject_names", colnames(x@subject_meta))
    unknown_vars <- setdiff(all.vars(subset_call), valid_vars)
    if (length(unknown_vars) > 0) {
        stop(sprintf(
            "subset references unknown variable(s): %s. Variables must be 'subject_names' or a column of subject_meta: %s",
            toString(unknown_vars), toString(colnames(x@subject_meta))
        ))
    }
    # Evaluate the subset expression in an environment of subject_meta columns plus subject_names
    eval_env <- list2env(c(as.list(x@subject_meta), list(subject_names = x@subject_names)))
    keep_subjects <- eval(subset_call, envir = eval_env)
    if (sum(keep_subjects, na.rm = TRUE) == 0) {
        warning("No subjects matched the subset criteria; returning an empty linkage_summary object.")
    }
    # subset the linkage_summary object, dropping unused factor levels so the
    # remaining subject_names factor accurately reflects the subjects kept
    new_subject_names <- droplevels(x@subject_names[keep_subjects])
    new_subject_meta <- x@subject_meta[keep_subjects, , drop = FALSE]
    new_subject_linkages <- x@subject_linkages[keep_subjects]
    new_link_sum <- linkage_summary(
        subject_names = new_subject_names,
        subject_meta = new_subject_meta,
        subject_linkages = new_subject_linkages
    )
    return(new_link_sum)
})
