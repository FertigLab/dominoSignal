#' Plot differential linkages among domino results ranked by a comparative statistic
#'
#' Plot differential linkages among domino results ranked by a comparative statistic
#'
#' @param differential_linkages a data frame output from the [test_differential_linkages()] function
#' @param test_statistic column name of differential_linkages where the test statistic used for ranking linkages is
#'   stored (ex. 'p.value')
#' @param stat_range a two value vector of the minimum and maximum values of test_statistic for
#'  plotting linkage features
#' @param stat_ranking 'ascending' (lowest value of test statisic is colored red and plotted at the top) or
#'   'descending' (highest value of test statistic is colored red and plotted at the top).
#' @param group_palette a named vector of colors to use for each group being compared
#' @return A heatmap-class object of features ranked by test_statistic annotated with the proportion of subjects
#'   that showed active linkage of the features.
#' @export
#' @examples
#' example(build_domino, echo = FALSE)
#' example(test_differential_linkages, echo = FALSE)
#' plot_differential_linkages(
#'     differential_linkages = tiny_differential_linkage_c1,
#'     test_statistic = "p.value",
#'     stat_ranking = "ascending"
#' )
#'
plot_differential_linkages <- function(
    differential_linkages, test_statistic, stat_range = c(0, 1),
    stat_ranking = c("ascending", "descending"), group_palette = NULL) {

    if (!test_statistic %in% colnames(differential_linkages)) {
        stop("test statistic '", test_statistic, "' not present in colnames(differential_linkages)")
    }
    if (identical(stat_ranking, c("ascending", "descending"))) {
        warning("stat_ranking order not specified. Defaulting to ascending order")
        stat_ranking <- "ascending"
    }
    if (!stat_ranking %in% c("ascending", "descending")) {
        stop("stat_ranking must be 'ascending' or 'descending'")
    }
    # limit to features within stat range
    dframe <- differential_linkages[differential_linkages[[test_statistic]] >= stat_range[1] &
            differential_linkages[[test_statistic]] <= stat_range[2], ]
    if (nrow(dframe) == 0) {
        stop("No features with '", test_statistic, "' within stat_range")
    }
    # order df by plot statistic
    if (stat_ranking == "ascending") {
        dframe <- dframe[order(dframe[[test_statistic]], dframe[["total_count"]], decreasing = FALSE), ]
        stat_gradient <- c("#FF0000", "#FFFFFF")
    }
    if (stat_ranking == "descending") {
        dframe <- dframe[order(dframe[[test_statistic]], dframe[["total_count"]], decreasing = TRUE), ]
        stat_gradient <- c("#FFFFFF", "#FF0000")
    }
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
