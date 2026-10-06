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

test_that("dom_info returns recorded package versions, or NULL when absent", {
    info <- dom_info(tiny_dom1)
    expect_named(info, c("create", "build", "create_variables", "build_variables", "create_version", "build_version"))
    expect_null(info$create_version)
    expect_null(info$build_version)

    dom <- tiny_dom1
    dom@misc$create_version <- "1.7.0"
    dom@misc$build_version <- "1.8.0"
    info <- dom_info(dom)
    expect_equal(info$create_version, "1.7.0")
    expect_equal(info$build_version, "1.8.0")
})

test_that("dom_info returns creation variables, or NULL when absent", {
    expect_null(dom_info(tiny_dom1)$create_variables)
    dom <- tiny_dom1
    dom@misc$create_vars <- list(tf_selection_method = "clusters")
    expect_equal(dom_info(dom)$create_variables, list(tf_selection_method = "clusters"))
})
