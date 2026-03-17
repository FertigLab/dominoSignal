test_that("test_differential_linkages runs", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- test_differential_linkages(
        linkage_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        test_name = "fishers.exact"
    )

    expect_s3_class(out, "data.frame")
})
