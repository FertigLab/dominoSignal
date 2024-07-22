test_that("Creating a circos plot data frame", {
  arc_df <- circos_lr_shape_data(
    dom = v0.2.1$pbmc_dom_built_tiny, 
    receptor="CXCR3", 
    ligand_expression_threshold = 0.01, cell_idents = NULL
  )
  expect_true(is(arc_df, "data.frame"))
  
})