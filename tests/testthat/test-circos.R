test_that("domino object interpretation by obtain_circos_expression", {
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
  
})