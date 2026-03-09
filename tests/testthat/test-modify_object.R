test_that("add_rl_column function runs", {

    

    map <- data.frame(a = c("x", "y"), b = c("A", "B"), stringsAsFactors = FALSE)
    conv <- data.frame(old = c("x", "y"), new = c("x1", "y1"), stringsAsFactors = FALSE)
    expect_no_error(add_rl_column(map = map, map_ref = "a", conv = conv, new_name = "mapped"))
})

test_that("rename_clusters function works correctly", {
    
    # Define the cluster conversion, Z/new_clust does not match data intentionally
    clust_conv <- c("W", "X", "Y", "Z")
    names(clust_conv) <- c(levels(tiny_dom1@clusters[1:3]), "new_clust")
    # Make sure it runs, to start
    expect_no_error(rename_clusters(tiny_dom1, clust_conv))

    # Run the function
    dom_renamed <- rename_clusters(tiny_dom1, clust_conv)

    # Check that the clusters were renamed correctly (no Zs, no new_clust)
    expect_equal(levels(dom_renamed@clusters), c("W", "X", "Y"))
    expect_equal(colnames(dom_renamed@clust_de), c("W", "X", "Y"))
    expect_named(dom_renamed@linkages$clust_tf, c("W", "X", "Y"))
    expect_named(dom_renamed@linkages$clust_rec, c("W", "X", "Y"))
    expect_named(dom_renamed@linkages$clust_incoming_lig, c("W", "X", "Y"))
    expect_named(dom_renamed@linkages$clust_tf_rec, c("W", "X", "Y"))
    expect_equal(colnames(dom_renamed@signaling), paste0("L_", c("W", "X", "Y")))
    expect_equal(rownames(dom_renamed@signaling), paste0("R_", c("W", "X", "Y")))
    expect_named(dom_renamed@cl_signaling_matrices, c("W", "X", "Y"))
    expect_equal(
        colnames(dom_renamed@cl_signaling_matrices[["W"]]),
        paste0("L_", c("W", "X", "Y"))
    )
    expect_equal(
        colnames(dom_renamed@cl_signaling_matrices[["X"]]),
        paste0("L_", c("W", "X", "Y"))
    )
    expect_equal(
        colnames(dom_renamed@cl_signaling_matrices[["Y"]]),
        paste0("L_", c("W", "X", "Y"))
    )
})
