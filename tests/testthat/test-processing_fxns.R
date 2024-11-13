test_that("build domino does not fail with no TFs with p-value below thd", {
   expect_no_error(pbmc_dom_built <- build_domino(
        dom = v0.2.1$pbmc_dom_tiny,
        min_tf_pval = 1e-20,
        max_tf_per_clust = Inf,
        max_rec_per_tf = Inf,
        rec_tf_cor_threshold = 1e-20,
        min_rec_percentage = Inf
  ))
})