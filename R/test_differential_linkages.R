#' @import methods
#'
NULL

#' Statistical test for differential linkages across multiple domino results
#'
#' Statistical test for differential linkages across multiple domino results
#'
#' @param linkage_summary a [linkage_summary()] object
#' @param cluster the name of the cell cluster being compared across multiple domino results
#' @param group.by the name of the column in `linkage_summary@subject_meta` by which to group subjects for counting.
#' @param linkage a stored linkage from the domino object.
#'   Can compare any of 'tfs', 'rec', 'incoming_lig', 'tfs_rec', or 'rec_lig'
#' @param subject_names a vector of subject_names from the linkage_summary to be compared.
#'   If NULL, all subject_names in the linkage summary are included in counting.
#' @param test_name the statistical test used for comparison.
#' \itemize{
#'  \item{'fishers.exact'} : Fisher's exact test for the dependence of the proportion of subjects with an
#'   active linkage in the cluster on which group the subject belongs to in the group.by variable.
#'   Provides an odds ratio, p-value, and a Benjamini-Hochberg FDR-adjusted p-value (p.adj) for each linkage tested.
#' }
#' @return A data frame of results from the test of the differential linkages. Rows correspond to each
#'   linkage tested. Columns correspond to:
#' \itemize{
#'  \item{'cluster'} : the name of the cell cluster being compared
#'  \item{'linkage'} : the type of linkage being compared
#'  \item{'group.by'} : the grouping variable
#'  \item{'test_name'} : the test used for comparison
#'  \item{'feature'} : individual linkages compared
#'  \item{'test statistics'} : test statistics provided are based on test method.
#'   'fishers.exact' provides a odds ratio, p-value, and fdr-adjusted p-value.
#'  \item{'total_count'} : total number of subjects where the linkage is active
#'  \item{'X_count'} : number of subjects in each category of group.by (X) where the linkage is active
#'  \item{'total_n'} : number of total subjects compared
#'  \item{'X_n'} : total number of subjects in each category of group.by (X)
#' }
#' @export
#' @examples
#' tiny_differential_linkage_c1 <- test_differential_linkages(
#'     linkage_summary = mock_linkage_summary(), cluster = "C1", group.by = "group",
#'     linkage = "rec", test_name = "fishers.exact"
#' )
#'
test_differential_linkages <- function(
    linkage_summary, cluster, group.by, linkage = "rec_lig", subject_names = NULL,
    test_name = "fishers.exact"
) {

    check_arg(linkage_summary, allow_class = "linkage_summary", allow_len = 1)
    check_arg(cluster, allow_class = "character", allow_len = 1)
    check_arg(group.by, allow_class = "character", allow_len = 1, allow_values = colnames(linkage_summary@subject_meta))
    check_arg(linkage, allow_class = "character", allow_len = 1,
        allow_values = c("tfs", "rec", "incoming_lig", "tfs_rec", "rec_lig"))
    check_arg(subject_names, allow_class = c("factor", "character", "NULL"))
    check_arg(test_name, allow_class = "character", allow_len = 1, allow_values = "fishers.exact")

    if (is.null(subject_names)) {
        subject_names <- linkage_summary@subject_names
    }
    # count the number of groups
    subject_count <- as.data.frame(table(linkage_summary@subject_meta[[group.by]]))
    colnames(subject_count) <- c(group.by, "total")
    group_levels <- subject_count[[group.by]]
    count_link <- count_linkage(
        linkage_summary = linkage_summary, cluster = cluster, linkage = linkage,
        group.by = group.by, subject_names = subject_names
    )
    # initiate data frame for storing results
    n <- nrow(count_link)
    result_df <- data.frame(cluster = rep(cluster, n), linkage = rep(linkage, n), group.by = rep(
        group.by,
        n
    ), test_name = rep(test_name, n), feature = count_link[["feature"]])
    # empty contigency table
    test_mat <- matrix(data = NA, nrow = nrow(subject_count), ncol = 2)
    rownames(test_mat) <- subject_count[[group.by]]
    colnames(test_mat) <- c("linkage_present", "linkage_absent")
    test_template <- as.data.frame(test_mat)
    if (test_name == "fishers.exact") {
        test_result <- as.data.frame(t(vapply(
            result_df[["feature"]],
            FUN.VALUE = numeric(2), FUN = function(x) {
                feat_count <- count_link[count_link[["feature"]] == x, !colnames(count_link) %in% c(
                    "feature",
                    "total_count"
                )]
                feat_count <- vapply(feat_count, FUN.VALUE = numeric(1), FUN = as.numeric)
                # fill contingency table
                test_df <- test_template
                test_df[["linkage_present"]] <- feat_count
                test_df[["linkage_absent"]] <- subject_count[["total"]] - feat_count
                # conduct test
                test <- fisher.test(test_df)
                odds.ratio <- test$estimate
                p.value <- test$p.value
                res <- c(odds.ratio, p.value)
                res <- setNames(res, c("odds.ratio", "p.value"))
                return(res)
            }
        )))
        # include fdr-adjusted p-values
        test_result[["p.adj"]] <- p.adjust(p = test_result[["p.value"]], method = "fdr")
    }
    # append test result columns to result
    result_df <- cbind(result_df, test_result)
    # append counts of active linkages for each group
    count_append <- count_link[, colnames(count_link) != "feature"]
    colnames(count_append) <- vapply(
        colnames(count_append),
        FUN.VALUE = character(1), FUN = function(x) {
            if (grepl("_count$", x)) {
                return(x)
            } else {
                return(paste0(x, "_count"))
            }
        }
    )
    result_df <- cbind(result_df, count_append)
    # append total counts of subjects in each group and the total number of subjects
    total_n <- sum(subject_count[["total"]])
    result_df[["total_n"]] <- total_n
    for (g in group_levels) {
        g_n <- as.numeric(subject_count[subject_count[[group.by]] == g, "total"])
        result_df[[paste0(g, "_n")]] <- g_n
    }
    return(result_df)
}
