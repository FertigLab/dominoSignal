test_that("build_domino runs with tiny object", {
    expect_s4_class(
        build_domino(
            dom = tiny_created_dom1,
            max_tf_pval = 0.05,
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
        max_tf_pval = 0,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 1e-20,
        min_rec_percentage = 1
    ))

    # TODO: Return warning/message to user if this happens:
})

# Check that build_domino runs correctly with custom rl_map with only gene and type columns
test_that("build_domino runs with custom rl_map with minimum required columns", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = TRUE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    dom_custom <- build_domino(dom_custom,
        max_tf_pval = 0.05,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    expect_false(all(dom_signaling(dom_custom) == 0))
})

# test that build_domino runs correctly with a custom rl_map with only gene and type
# columns when use_complexes = FALSE
test_that("build_domino runs with custom rl_map with minimum required columns and use_complexes = FALSE", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = FALSE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    dom_custom <- build_domino(dom_custom,
        max_tf_pval = 0.05,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    expect_false(all(dom_signaling(dom_custom) == 0))
})

# When only one receptor passes expression threshold, build_domino should still identify it.
test_that("build_domino identifies an expressed receptor when only one receptor passes rec_min_thresh", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = FALSE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    expect_equal(nrow(dom_custom@misc$cl_rec_percent), 1)
    expect_equal(rownames(dom_custom@misc$cl_rec_percent), "CXCR3")

    dom_custom <- build_domino(dom_custom,
        max_tf_pval = 0.05,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    clust_rec <- dom_linkages(dom_custom, "receptor", by_cluster = TRUE)
    expect_true("CXCR3" %in% clust_rec[["CD8_T_cell"]])
    expect_true("CXCR3" %in% clust_rec[["B_cell"]])
    expect_false("CXCR3" %in% clust_rec[["CD14_monocyte"]])
})

# When cluster receives only one incoming ligand, build_domino should still compute nonzero signaling
test_that("build_domino computes non-zero signaling for a cluster with a single incoming ligand", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = FALSE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    dom_custom <- build_domino(dom_custom,
        max_tf_pval = 0.05,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    incoming_ligs <- dom_linkages(dom_custom, "incoming-ligand", by_cluster = TRUE)
    expect_length(incoming_ligs[["CD8_T_cell"]], 1)

    expected <- vapply(levels(tiny_clusters1), function(cl) {
        max(mean(tiny_zscores1["CCL20", tiny_clusters1 == cl]), 0)
    }, numeric(1))
    observed <- unname(unlist(dom_signaling(dom_custom)["R_CD8_T_cell", ]))
    expect_equal(observed, unname(expected), tolerance = 1e-8)
})

# When complexes are used without name columns, signaling should still be produced
test_that("build_domino runs with rl_map without name columns and use_complexes = TRUE", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = TRUE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    expect_named(dom_custom@linkages$complexes, c("ITGB4,ITGA6", "IL7R,IL2RG"))

    dom_custom <- build_domino(dom_custom,
        max_tf_pval = 0.05,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    expect_false(all(dom_signaling(dom_custom) == 0))
})

test_that("build_domino returns no active TFs when no p-value passes max_tf_pval", {
    expect_message(built <- build_domino(
        dom = tiny_created_dom1, max_tf_pval = 0, max_tf_per_clust = Inf, max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01
    ), "No TFs in cluster .* will be called active as no p-values pass the max_tf_pval threshold.*")
    expect_true(all(lengths(built@linkages$clust_tf) == 0))
    expect_true(all(lengths(built@linkages$clust_rec) == 0))
    expect_true(all(built@signaling == 0))
})

test_that("build_domino averages complex ligand components when a cluster has a single incoming complex ligand", {
    built <- build_domino(
        dom = tiny_created_dom1, max_tf_pval = 0.05, max_tf_per_clust = Inf, max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01
    )
    expect_identical(built@linkages$clust_incoming_lig[["CD14_monocyte"]], "integrin_a6b4_complex")
    mono_mat <- built@cl_signaling_matrices[["CD14_monocyte"]]
    expect_identical(rownames(mono_mat), "integrin_a6b4_complex")

    # floored cluster mean z-score of each component gene, averaged across components
    comp_means <- vapply(levels(tiny_clusters1), function(cl) {
        mean(pmax(rowMeans(tiny_zscores1[c("ITGB4", "ITGA6"), tiny_clusters1 == cl, drop = FALSE]), 0))
    }, numeric(1))
    expect_equal(unname(mono_mat[1, ]), unname(comp_means))
    expect_equal(unname(built@signaling["R_CD14_monocyte", ]), unname(comp_means))
})
