test_that("plot_network functions run", {

    expect_no_error(signaling_network(tiny_dom1, edge_weight = 0.1, scale = "none", normalize = "none"))
    expect_no_error(gene_network(tiny_dom1, clust = levels(tiny_dom1@clusters)[1], layout = "grid"))
})
