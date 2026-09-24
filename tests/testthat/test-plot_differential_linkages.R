# a synthetic differential linkage result with a controlled range of test statistics and one
# feature (F5) where odds.ratio is NA, mimicking a mix of 2-group and >2-group comparisons
mock_diff_linkage <- data.frame(
    cluster = rep("C1", 5),
    linkage = rep("rec", 5),
    group.by = rep("group", 5),
    test_name = rep("fishers.exact", 5),
    feature = c("F1", "F2", "F3", "F4", "F5"),
    odds.ratio = c(0.5, 1.5, 2.5, 3.5, NA),
    p.value = c(0.01, 0.2, 0.03, 0.5, 0.8),
    p.adj = c(0.02, 0.3, 0.04, 0.6, 0.9),
    total_count = c(4, 3, 5, 2, 1),
    A_count = c(2, 1, 3, 1, 1),
    B_count = c(2, 2, 2, 1, 0),
    total_n = rep(10, 5),
    A_n = rep(5, 5),
    B_n = rep(5, 5),
    stringsAsFactors = FALSE
)

# get the continuous color mapping function used to color heatmap cells
get_col_fun <- function(plt) {
    plt@ht_list[[1]]@matrix_color_mapping@col_fun
}

# convert a color to the hex+alpha format returned by circlize::colorRamp2
expected_hex <- function(col) {
    rgb_vals <- grDevices::col2rgb(col)[, 1]
    sprintf("#%02X%02X%02XFF", rgb_vals[1], rgb_vals[2], rgb_vals[3])
}

test_that("plot_differential_linkages runs", {
    plt <- plot_differential_linkages(
        differential_linkages = tiny_differential_linkage,
        test_statistic = "p.value",
        stat_ranking = "ascending"
    )

    expect_s4_class(plt, "HeatmapList")
})

test_that("plot_differential_linkages orders features by stat_ranking", {
    plt_asc <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "p.value",
        stat_ranking = "ascending"
    )
    expect_identical(rownames(plt_asc@ht_list[[1]]@matrix), c("F1", "F3", "F2", "F4", "F5"))

    plt_desc <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "p.value",
        stat_ranking = "descending"
    )
    expect_identical(rownames(plt_desc@ht_list[[1]]@matrix), c("F5", "F4", "F2", "F3", "F1"))
})

test_that("plot_differential_linkages filters features by stat_range", {
    plt <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "odds.ratio",
        stat_range = c(1, 3),
        stat_ranking = "ascending"
    )
    expect_setequal(rownames(plt@ht_list[[1]]@matrix), c("F2", "F3"))

    expect_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "odds.ratio",
            stat_range = c(10, 20)
        ),
        "No features with 'odds.ratio' within stat_range"
    )
})

test_that("stat_range is not restricted to c(0, 1) for non p-value statistics", {
    expect_no_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "odds.ratio",
            stat_range = c(0, 4)
        )
    )
    # p.value and p.adj are still restricted to the c(0, 1) range
    expect_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "p.value",
            stat_range = c(-1, 2)
        )
    )
})

test_that("rows with NA test_statistic values are dropped rather than plotted as blank rows", {
    plt <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "odds.ratio",
        stat_range = c(0, 4)
    )
    mat <- plt@ht_list[[1]]@matrix
    expect_false("F5" %in% rownames(mat))
    expect_equal(nrow(mat), 4)
    expect_false(anyNA(mat))
})

test_that("an error is raised when all test_statistic values are NA", {
    all_na_linkage <- mock_diff_linkage
    all_na_linkage$odds.ratio <- NA_real_

    expect_error(
        plot_differential_linkages(
            differential_linkages = all_na_linkage,
            test_statistic = "odds.ratio"
        ),
        "All values in test_statistic column 'odds.ratio' are NA; cannot plot."
    )
})

test_that("a missing test_statistic column gives an informative error, not the NA message", {
    expect_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "not_a_column"
        ),
        "Required variables not_a_column not found"
    )
})

test_that("the default gradient colors the minimum and maximum of stat_range", {
    plt <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "p.value",
        stat_range = c(0, 1)
    )
    col_fun <- get_col_fun(plt)
    expect_equal(col_fun(0), expected_hex("red"))
    expect_equal(col_fun(1), expected_hex("gray90"))
})

test_that("a custom gradient is applied regardless of the order names are given in", {
    plt <- plot_differential_linkages(
        differential_linkages = mock_diff_linkage,
        test_statistic = "odds.ratio",
        stat_range = c(0, 4),
        gradient = c(maximum = "blue", minimum = "yellow")
    )
    col_fun <- get_col_fun(plt)
    expect_equal(col_fun(0), expected_hex("yellow"))
    expect_equal(col_fun(4), expected_hex("blue"))
})

test_that("gradient must be a named vector with 'minimum' and 'maximum' names", {
    expect_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "p.value",
            gradient = c("red", "gray90")
        ),
        "names"
    )
    expect_error(
        plot_differential_linkages(
            differential_linkages = mock_diff_linkage,
            test_statistic = "p.value",
            gradient = c(low = "red", high = "gray90")
        ),
        "gradient must be a named vector with names 'minimum' and 'maximum'"
    )
})
