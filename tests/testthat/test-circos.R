test_that("domino object interpretation by obtain_circos_expression", {
  # start temporary graphics device for testing to preserve package enviroment
  png(filename = paste0(tempdir(), "/", "ts.png"))
  
  # testing domino object
  data(CellPhoneDB)
  rl_map_tiny <- create_rl_map_cellphonedb(
    genes = CellPhoneDB$genes_tiny,
    proteins = CellPhoneDB$proteins_tiny,
    interactions = CellPhoneDB$interactions_tiny,
    complexes = CellPhoneDB$complexes_tiny)
  
  data(SCENIC)
  regulon_list_tiny <- create_regulon_list_scenic(
    regulons = SCENIC$regulons_tiny)
  
  data(PBMC)
  pbmc_dom_tiny <- create_domino(
    rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
    counts = PBMC$RNA_count_tiny, z_scores = PBMC$RNA_zscore_tiny,
    clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
    use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE)
  
  dom <- build_domino(
    dom = pbmc_dom_tiny, min_tf_pval = .05, max_tf_per_clust = Inf,
    max_rec_per_tf = Inf, rec_tf_cor_threshold = .1, min_rec_percentage = 0.01
  )
  
  # parent function runs to completion without error
  expect_no_error(circos_ligand_receptor(dom, receptor = "CXCR3"))
  
  # specification of a single ligand
  circos_df_ccl20 <- obtain_circos_expression(dom, receptor = "CXCR3", ligands = "CCL20")
  expect_equal(unique(circos_df_ccl20$ligand), "CCL20")
  
  # fail without specification of a ligands
  expect_error(obtain_circos_expression(dom, receptor = "CXCR3"))
  
  dev.off()
})

test_that("Cell types with hyphenated names can be plotted", {
  # start temporary graphics device for testing to preserve package enviroment
  png(filename = paste0(tempdir(), "/", "ts.png"))
  
  ts_signaling_df <- data.frame(
    origin = c("CT_1-L1", "CT_1-L2", "CT_2-L1", "CT_2-L2", "CT_3-L1", "CT_3-L2"),
    destination = "R1",
    mean.expression = c(0,0.5, 1,0.5, 0.5,0.5),
    sender = c("CT_1", "CT_1", "CT_2", "CT_2", "CT_3", "CT_3"),
    ligand = c("L1", "L2", "L1", "L2", "L1", "L2"),
    receptor = "R1",
    scaled.mean.expression = c(0,0.5, 1,0.5, 0.5,0.5),
    ligand.arc = 1,
    receptor.arc = 4/6
  )
  
  # plot can render from a manually written signaling data.frame
  expect_no_error(render_circos_ligand_receptor(ts_signaling_df, receptor = "R1"))
  
  # plot can render if all cell type underscores are replaced by hyphens
  ts_signaling_df_hyphen <- ts_signaling_df
  ts_signaling_df_hyphen$origin <- gsub("_", "-", ts_signaling_df$origin)
  ts_signaling_df_hyphen$sender <- gsub("_", "-", ts_signaling_df$sender)
  
  expect_no_error(render_circos_ligand_receptor(ts_signaling_df_hyphen, receptor = "R1"))
  
  dev.off()
})
