test_that("count_linkage runs", {

    linkage_sum <- tiny_linkage_summary
    clust <- names(linkage_sum@subject_linkages[[1]])[1]

    out <- count_linkage(
        linkage_summary = linkage_sum,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )

    expect_s3_class(out, "data.frame")
})
