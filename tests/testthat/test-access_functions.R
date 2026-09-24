test_that("access functions run", {
    expect_no_error(dom_database(tiny_dom1))
    expect_no_error(dom_zscores(tiny_dom1))
    expect_no_error(dom_counts(tiny_dom1))
    expect_no_error(dom_clusters(tiny_dom1))
    expect_no_error(dom_tf_activation(tiny_dom1))
    expect_no_error(dom_correlations(tiny_dom1))
    expect_no_error(dom_linkages(tiny_dom1, "receptor-ligand"))
    expect_no_error(dom_signaling(tiny_dom1))
    expect_no_error(dom_de(tiny_dom1))
    expect_no_error(dom_info(tiny_dom1))
    expect_no_error(dom_network_items(tiny_dom1))
})

test_that("dom_signalling is a synonym of dom_signaling", {
    expect_identical(dom_signalling, dom_signaling)
    expect_identical(dom_signalling(tiny_dom1), dom_signaling(tiny_dom1))
    expect_identical(
        dom_signalling(tiny_dom1, cluster = "CD14_monocyte"),
        dom_signaling(tiny_dom1, cluster = "CD14_monocyte")
    )
})
