test_that("test_differential_linkages runs", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- test_differential_linkages(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )

    expect_s3_class(out, "data.frame")
})
