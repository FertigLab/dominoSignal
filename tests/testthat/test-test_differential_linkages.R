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

test_that("test_differential_linkages defaults to include all subjects", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    out_full <- test_differential_linkages(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec"
    )
    out_restricted <- test_differential_linkages(
        link_summary = tiny_linkage_summary,
        cluster = clust,
        group.by = "group",
        linkage = "rec",
        subject_names = c("dom1", "dom2", "dom3")
    )

    expect_equal(out_full, out_restricted)
    expect_equal(out_full$total_n, c(3, 3))
    expect_equal(out_full$A_n, c(2, 2))
    expect_equal(out_full$B_n, c(1, 1))
})

test_that("test_differential_linkages counts reflect only the included subjects, not the full link_summary", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]
    tiny_linksum3 <- tiny_linkage_summary
    tiny_linksum3@subject_meta$group2 <- c("A", "B", "C")

    out <- test_differential_linkages(
        link_summary = tiny_linksum3,
        cluster = clust,
        group.by = "group2",
        linkage = "rec",
        subject_names = c("dom1", "dom2")
    )

    expect_equal(out$total_n, c(2, 2))
    expect_equal(out$A_n, c(1, 1))
    expect_equal(out$B_n, c(1, 1))
    expect_false("C_n" %in% colnames(out))
})

test_that("test_differential_linkages errors when subject_names leaves fewer than 2 groups", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    expect_error(
        test_differential_linkages(
            link_summary = tiny_linkage_summary,
            cluster = clust,
            group.by = "group",
            linkage = "rec",
            subject_names = c("dom1", "dom2")
        ),
        "At least 2 groups"
    )
})

test_that("test_differential_linkages errors on subject_names not present in link_summary", {
    clust <- names(tiny_linkage_summary@subject_linkages[[1]])[1]

    expect_error(
        test_differential_linkages(
            link_summary = tiny_linkage_summary,
            cluster = clust,
            group.by = "group",
            linkage = "rec",
            subject_names = c("dom1", "bogus")
        ),
        "bogus"
    )
})
