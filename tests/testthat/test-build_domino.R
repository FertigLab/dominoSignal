test_that("build_domino runs with tiny object", {
    expect_s4_class(
        build_domino(
            dom = tiny_created_dom1,
            min_tf_pval = 0.05,
            max_tf_per_clust = 3,
            max_rec_per_tf = 3,
            rec_tf_cor_threshold = 0.1,
            min_rec_percentage = 0.01
        ),
        "domino"
    )
})

test_that("build_domino does not fail with no TFs with p-value below threshold", {
    expect_no_error(build_domino(
        dom = tiny_created_dom1,
        min_tf_pval = 0,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 1e-20,
        min_rec_percentage = Inf
    ))

    # TODO: Return warning/message to user if this happens:
})