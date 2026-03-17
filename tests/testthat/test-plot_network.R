test_that("signaling_network function runs", {

    expect_no_error(signaling_network(tiny_dom1, edge_weight = 0.1, scale = "none", normalize = "none"))
})

test_that("signaling_network returns error if no signaling is found", {
    expect_error(signaling_network(tiny_dom1, showOutgoingSignalingClusts = "CD14_monocyte",
            showIncomingSignalingClusts = "CD8_T_cell"), regexp = "No signaling found")
})

test_that("signaling_network returns a graph with scaled vertices even if scaling results in NA values", {
    expect_no_error(signaling_network(tiny_dom1, showOutgoingSignalingClusts = "CD8_T_cell",
            scale_by = "lig_sig"))
    expect_no_error(signaling_network(tiny_dom1, showIncomingSignalingClusts = "B_cell"))
})

test_that("gene_network function runs", {
    expect_no_error(gene_network(tiny_dom1, clust = levels(tiny_dom1@clusters)[1], layout = "grid"))
})

test_that("gene_network can handle multiple clusters for incoming and outgoing signaling", {
    expect_no_error(gene_network(tiny_dom1, clust = levels(tiny_dom1@clusters)[1:2], 
            OutgoingSignalingClust = levels(tiny_dom1@clusters)[2:3], layout = "grid"))
    
    gn_plot <- gene_network(tiny_dom1, clust = levels(tiny_dom1@clusters)[1:2], 
        OutgoingSignalingClust = levels(tiny_dom1@clusters)[2:3], layout = "grid")
    expect_named(gn_plot, c("graph", "layout"))
    expect_s3_class(gn_plot[[1]], "igraph")
    expect_contains(class(gn_plot[[2]]), "matrix")
})
