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

test_that("signalling_network is a synonym of signaling_network", {
    expect_identical(signalling_network, signaling_network)
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

test_that("signaling_network handles normalization of clusters with no signaling", {
    dom <- tiny_dom1
    dom@signaling["R_B_cell", ] <- 0
    expect_no_error(signaling_network(dom, normalize = "rec_norm"))
    expect_no_error(signaling_network(dom, normalize = "lig_norm"))
})

test_that("signaling_network handles cluster names containing 'L_' or 'R_'", {
    dom <- rename_clusters(tiny_dom1, c(CD8_T_cell = "CD8_T_CELL_1", CD14_monocyte = "MONO_R_1", B_cell = "B"))
    # vertex names involved in colors and ligand sizing; only prefix should be stripped
    expect_no_warning(signaling_network(dom, scale_by = "lig_sig"))
    expect_no_error(signaling_network(dom, cols = c(CD8_T_CELL_1 = "red", MONO_R_1 = "blue", B = "green")))
})

test_that("gene_network requires clust and a built object", {
    expect_error(gene_network(tiny_dom1), "Please provide clust as one or more receptor clusters in the domino object")
    expect_error(gene_network(tiny_created_dom1, clust = "B_cell"), "Please build a signaling network")
})

test_that("gene_network sizes each ligand once when several receiving clusters share it", {
    clusts <- levels(tiny_dom1@clusters)[c(1, 3)]
    outgoing <- levels(tiny_dom1@clusters)
    net <- gene_network(tiny_dom1, clust = clusts, OutgoingSignalingClust = outgoing, lig_scale = 1)
    sizes <- stats::setNames(igraph::V(net$graph)$size, names(igraph::V(net$graph)))
    # ligand expression rows are identical across receiving clusters, so one matrix gives the expected size
    ligs <- intersect(names(sizes), rownames(tiny_dom1@cl_signaling_matrices[[clusts[1]]]))
    expect_gt(length(ligs), 0)
    mat <- tiny_dom1@cl_signaling_matrices[[clusts[1]]][ligs, paste0("L_", outgoing), drop = FALSE]
    expect_equal(sizes[ligs], 0.5 * rowSums(mat))
})