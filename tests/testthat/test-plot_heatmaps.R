test_that("plot_heatmaps functions run", {
    rec_clust <- levels(tiny_dom1@clusters)[1]

    expect_no_error(signaling_heatmap(tiny_dom1))
    expect_no_error(incoming_signaling_heatmap(tiny_dom1, rec_clust = rec_clust))
    expect_no_error(feat_heatmap(tiny_dom1, title = FALSE))
    expect_no_error(cor_heatmap(tiny_dom1, title = FALSE))
})

test_that("signalling_heatmap is a synonym of signaling_heatmap", {
    expect_identical(signalling_heatmap, signaling_heatmap)
})

test_that("incoming_signalling_heatmap is a synonym of incoming_signaling_heatmap", {
    expect_identical(incoming_signalling_heatmap, incoming_signaling_heatmap)
})
