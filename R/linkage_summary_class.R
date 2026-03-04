#' The domino linkage summary class
#'
#' The linkage summary class contains linkages established in multiple domino
#' objects through gene regulatory network inference and reference to receptor-
#' ligand data bases. A data frame summarizing meta features that describe the
#' domino objects compared in the linkage summary facilitates comparisons of
#' established linkages and differential signaling interactions across categorical
#' sample covariates.
#'
#' @slot subject_names unique names for each domino result included in the summary
#' @slot subject_meta data.frame with each row describing one subject and columns describing features of the subjects by which to draw comparisons of signaling networks
#' @slot subject_linkages nested list of linkages inferred for each subject. Lists are stored in a heirarchical structure of subject-cluster-linkage where linkages include transcription factors (tfs), linkages between transcription factors and receptors (tfs_rec), active receptors (rec), possible receptor-ligand interactions (rec_lig), and incoming ligands (incoming_lig)
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