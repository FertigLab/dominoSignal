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

test_that("test_differential_linkages runs with > 2 groups", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]
    tiny_linksum3 <- tiny_linkage_summary
    tiny_linksum3@subject_meta$group2 <- c("A", "B", "C")

    expect_no_error(out <- test_differential_linkages(
        link_summary = tiny_linksum3,
        cluster = clust,
        group.by = "group2",
        linkage = "rec"
    ))

    expect_s3_class(out, "data.frame")
    expect_all_equal(out$odds.ratio, NA_real_)
})
