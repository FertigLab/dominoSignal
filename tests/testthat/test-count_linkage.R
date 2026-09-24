test_that("count_linkage runs", {

    linkage_sum <- tiny_linkage_summary
    clust <- names(linkage_sum@subject_linkages[[1]])[1]

    out <- count_linkage(
        link_summary = linkage_sum,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )

    expect_s3_class(out, "data.frame")
})

test_that("count_linkage counts across all subjects by default", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )

    expect_equal(out$feature, c("CXCR3", "IL7_receptor"))
    expect_equal(out$total_count, c(3, 3))
    expect_equal(out$A, c(2, 2))
    expect_equal(out$B, c(1, 1))
})

test_that("count_linkage restricts counting to the given subject_names", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        subject_names = c("dom1", "dom2")
    )

    expect_equal(out$feature, c("CXCR3", "IL7_receptor"))
    expect_equal(out$total_count, c(2, 2))
    expect_equal(out$A, c(2, 2))
    expect_false("B" %in% colnames(out))
})

test_that("count_linkage restricts counting to a single subject_names with no group.by", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        linkage = "rec",
        subject_names = "dom1"
    )

    expect_equal(colnames(out), c("feature", "total_count"))
    expect_equal(out$total_count, c(1, 1))
})

test_that("count_linkage output is order-independent with respect to subject_names input order", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out_forward <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        subject_names = c("dom1", "dom2")
    )
    out_reversed <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        subject_names = c("dom2", "dom1")
    )

    expect_equal(out_forward, out_reversed)
})

test_that("count_linkage accepts subject_names as a factor", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        linkage = "rec",
        subject_names = factor(c("dom1", "dom2"))
    )

    expect_equal(out$total_count, c(2, 2))
})

test_that("count_linkage matches the default when subject_names names every subject", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out_default <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )
    out_explicit <- count_linkage(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        subject_names = c("dom1", "dom2", "dom3")
    )

    expect_equal(out_default, out_explicit)
})

test_that("count_linkage errors on subject_names not present in link_summary", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    expect_error(
        count_linkage(
            link_summary = tiny_linkage_summary,
            cluster = clust,
            linkage = "rec",
            subject_names = c("dom1", "bogus")
        ),
        "bogus"
    )
})
