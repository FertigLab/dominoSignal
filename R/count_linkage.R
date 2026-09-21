#' Count occurrences of linkages across multiple domino results from a linkage summary
#' 
#' Count occurrences of linkages across multiple domino results from a linkage summary
#' 
#' @param link_summary a [linkage_summary()] object
#' @param cluster the name of the cell cluster being compared across multiple domino results
#' @param group.by the name of the column in `linkage_summary@subject_meta` by which to group subjects for counting.
#'  If NULL, only total counts of linkages for linkages in the cluster across all subjects is given.
#' @param linkage a stored linkage from the domino object.
#'  Can compare any of 'tfs', 'rec', 'incoming_lig', 'tfs_rec', or 'rec_lig'
#' @param subject_names a vector of subject_names from the linkage_summary to restrict counting to.
#'   If NULL, all subjects in `link_summary` are included.
#' @return A data frame with columns for the unique linkage features and the counts of how many times the linkage
#'  occurred across the compared domino results. If group.by is used, counts of the linkages are also provided as
#'  columns named by the unique values of the group.by variable.
#' @export
#' @seealso [summarize_linkages()] to create [linkage_summary()] objects, [test_differential_linkages()] to compare links across groups
#' @family differentials
#' @examples
#' data(LinkageSummary)
#' count_linkage(
#'   link_summary = LinkageSummary$linkage_sum_tiny, cluster = "C1", 
#'   group.by = "group", linkage = "rec")
#' count_linkage(
#'   link_summary = LinkageSummary$linkage_sum_tiny, cluster = "C1", 
#'   group.by = "group", linkage = "tfs_rec", subject_names = c("P1", "P2", "P4", "P5"))

count_linkage <- function(link_summary, cluster, group.by = NULL, linkage = "rec_lig", subject_names = NULL) {
    check_arg(link_summary, allow_class = "linkage_summary", allow_len = 1)
    check_arg(cluster, allow_class = "character", allow_len = 1)
    check_arg(group.by, allow_class = c("character", "NULL"), allow_len = c(0, 1))
    check_arg(linkage, allow_class = "character", allow_len = 1,
        allow_values = c("tfs", "rec", "incoming_lig", "tfs_rec", "rec_lig"))
    check_arg(subject_names, allow_class = c("factor", "character", "NULL"))

    link_summary <- filter_by_subject_names(link_summary, subject_names)
    all_int_ls <- lapply(link_summary@subject_linkages, FUN = function(x) {
        return(x[[cluster]][[linkage]])
    })
    all_int <- unlist(all_int_ls)
    feature <- table(unlist(all_int))
    dframe <- data.frame(feature = names(feature), total_count = as.numeric(feature))
    if (!is.null(group.by)) {
        if (!group.by %in% colnames(link_summary@subject_meta)) {
            stop("group.by variable not present in subject_meta")
        }
        groups <- levels(factor(link_summary@subject_meta[[group.by]]))
        for (g in groups) {
            g_index <- link_summary@subject_meta[[group.by]] == g
            g_subjects <- link_summary@subject_meta[g_index, 1]
            int_count <- list()
            for (f in dframe[["feature"]]) {
                count <- vapply(g_subjects, FUN.VALUE = logical(1), FUN = function(x) {
                    f %in% link_summary@subject_linkages[[x]][[cluster]][[linkage]]
                })
                int_count[[f]] <- sum(count)
            }
            dframe[[g]] <- unlist(int_count)
        }
    }
    return(dframe)
}