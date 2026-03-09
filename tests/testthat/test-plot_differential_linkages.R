test_that("plot_differential_linkages runs", {
    plt <- plot_differential_linkages(
        differential_linkages = tiny_differential_linkage,
        test_statistic = "p.value",
        stat_ranking = "ascending"
    )

    expect_s4_class(plt, "HeatmapList")
})
