#' Plot differential linkages among domino results ranked by a comparative statistic
#'
#' Plot differential linkages among domino results ranked by a comparative statistic
#'
#' @param differential_linkages a data frame output from the [test_differential_linkages()] function
#' @param test_statistic column name of differential_linkages where the test statistic used for ranking linkages is
#'   stored (p.value, p.adj, odds.ratio)
#' @param stat_range a two value vector of the minimum and maximum values of test_statistic for
#'  plotting linkage features
#' @param stat_ranking whether to rank the features by ascending or descending values of test_statistic
#' @param group_palette a named vector of colors to use for each group being compared
#' @param gradient a named two value vector of colors to use for the minimum and maximum values of test_statistic
#' @return A heatmap-class object of features ranked by test_statistic annotated with the proportion of subjects
#'   that showed active linkage of the features.
#' @export
#' @family plotting
#' @family differentials
#' @examples
#' data(LinkageSummary)
#' plot_differential_linkages(
#'     differential_linkages = LinkageSummary$linkage_diff_tiny,
#'     test_statistic = "p.value",
#'     stat_range = c(0, 1),
#'     stat_ranking = "ascending"
#' )
#'
plot_differential_linkages <- function(
    differential_linkages, test_statistic, stat_range = c(0, 1),
    stat_ranking = c("ascending", "descending"), group_palette = NULL, gradient = c(minimum = "red", maximum = "gray90")) {
    
    stat_ranking <- match.arg(stat_ranking)

    check_arg(test_statistic, allow_class = "character", allow_len = 1)
    check_arg(differential_linkages, allow_class = "data.frame", need_vars = test_statistic)
    if (all(is.na(differential_linkages[[test_statistic]]))) {
        stop("All values in test_statistic column '", test_statistic, "' are NA; cannot plot.")
    }
    if (test_statistic %in% c("p.value", "p.adj")) {
        check_arg(stat_range, allow_class = "numeric", allow_len = 2, allow_range = c(0, 1))
    } else {
        check_arg(stat_range, allow_class = "numeric", allow_len = 2)
    }
    check_arg(gradient, allow_class = "character", allow_len = 2, need_names = TRUE)
    if (!setequal(names(gradient), c("minimum", "maximum"))) {
        stop("gradient must be a named vector with names 'minimum' and 'maximum'")
    }

    # limit to features within stat range, dropping any rows where test_statistic is NA
    stat_col <- differential_linkages[[test_statistic]]
    dframe <- differential_linkages[!is.na(stat_col) & stat_col >= stat_range[1] &
            stat_col <= stat_range[2], ]
    if (nrow(dframe) == 0) {
        stop("No features with '", test_statistic, "' within stat_range")
    }
    # order df by plot statistic
    if (stat_ranking == "ascending") {
        dframe <- dframe[order(dframe[[test_statistic]], dframe[["total_count"]], decreasing = FALSE), ]
    }
    if (stat_ranking == "descending") {
        dframe <- dframe[order(dframe[[test_statistic]], dframe[["total_count"]], decreasing = TRUE), ]
    }
    stat_gradient <- c(gradient["minimum"], gradient["maximum"])
    # values from test result for plotting
    cluster <- unique(dframe[["cluster"]])
    g_names_full <- colnames(dframe)[grepl("_n$", colnames(dframe)) & !grepl("^total_", colnames(dframe))]
    g_names <- gsub("_n", "", g_names_full, fixed = TRUE)
    # proportion bar for linkage feature in all subjects
    ha_subject <- ComplexHeatmap::HeatmapAnnotation(
        subjects = ComplexHeatmap::anno_barplot(matrix(ncol = 2, c(
            dframe[["total_count"]],
            dframe[["total_n"]] - dframe[["total_count"]]
        )), gp = grid::gpar(fill = c("black", "white"))), which = "row",
        annotation_name_gp = grid::gpar(fontsize = 8)
    )
    ha_subject@anno_list$subjects@label <- "All\nSubjects"
    # row annotation of linkage feature names
    ha_name <- ComplexHeatmap::rowAnnotation(feat = ComplexHeatmap::anno_text(dframe[["feature"]],
            location = 0, rot = 0))
    # plotted statistic for ordering results
    mat <- matrix(dframe[[test_statistic]], ncol = 1)
    rownames(mat) <- dframe[["feature"]]
    plotted <- ComplexHeatmap::Heatmap(matrix = mat, cluster_rows = FALSE, left_annotation = ha_name,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = grid::gpar(fontsize = 6))
        },
        column_title = paste0(cluster, ": ", test_statistic), name = test_statistic,
        col = circlize::colorRamp2(breaks = stat_range, colors = stat_gradient),
        height = nrow(mat) * grid::unit(0.25, "in"), width = grid::unit(1, "in")) + ha_subject
    # generate an heatmap annotation for each category
    if (is.null(group_palette)) {
        group_palette <- ggplot_col_gen(length(g_names))
        names(group_palette) <- g_names
    }
    for (i in seq_along(g_names)) {
        g <- g_names[i]
        g_count <- paste0(g, "_count")
        g_n <- paste0(g, "_n")
        ha <- ComplexHeatmap::HeatmapAnnotation(
            group = ComplexHeatmap::anno_barplot(matrix(ncol = 2, c(dframe[[g_count]],
                        dframe[[g_n]] - dframe[[g_count]])),
                gp = grid::gpar(fill = c(group_palette[g], "#FFFFFF"))),
            name = g, which = "row", annotation_name_gp = grid::gpar(fontsize = 8)
        )
        ha@anno_list$group@label <- g
        plotted <- plotted + ha
    }
    return(plotted)
}
