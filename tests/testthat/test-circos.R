test_that("Creating a circos plot data frame", {
  arc_df <- circos_lr_shape_data(
    dom = v0.2.1$pbmc_dom_built_tiny, 
    receptor = "CXCR3", ligands = "CCL20",
    ligand_expression_threshold = 0.01, cell_idents = NULL
  )
  expect_true(is(arc_df, "data.frame"))
  
  # subsetting to valid clusters does not cause a warning
  expect_no_warning(
    {
      # vector of multiple valid cell type
      circos_lr_shape_data(
        dom = v0.2.1$pbmc_dom_built_tiny, 
        receptor="CXCR3", 
        cell_idents = c("CD14_monocyte","CD3_T_cell")
      )
      # single character value of one valid cell type
      circos_lr_shape_data(
        dom = v0.2.1$pbmc_dom_built_tiny, 
        receptor="CXCR3", 
        cell_idents = "CD14_monocyte"
      )
    }
  )
  
  # providing invalid cluster names generates a warning
  expect_warning(
    circos_lr_shape_data(
      dom = v0.2.1$pbmc_dom_built_tiny, 
      receptor="CXCR3", 
      cell_idents = c("foobar")
    )
  )
})

# arc_df_2 <- circos_lr_shape_data(
#   dom = dominoSignal:::v0.2.1$pbmc_dom_built_tiny,
#   ligands = "CCL20",
#   receptor = "CXCR3", 
#   ligand_expression_threshold = 0.01, cell_idents = NULL
# )
# 
# circos_lr_plot(
#   arc_df = arc_df_2, ligands = "CCL20", receptor = "CXCR3"
# )
# 
# new_circos_ligand_receptor(
#   dom = dominoSignal:::v0.2.1$pbmc_dom_built_tiny, receptor = "CXCR3",
#   cell_idents = c("B_cell", "CD14_monocyte")
# )
# 
# pdf(file = "test.pdf", width = 8, height = 8)
# new_circos_ligand_receptor(
#   dom = dominoSignal:::v0.2.1$pbmc_dom_built_tiny, receptor = "CXCR3",
#   cell_idents = c("CD14_monocyte", "B_cell")
# )
# dev.off()

