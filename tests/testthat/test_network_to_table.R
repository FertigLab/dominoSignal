# Write tests for functions and component functions turning network
# data into a data frame.
# Consider adding results to sysdata for clearer code

library(testthat)

# Start with resolve_names
test_that("resolve_names replaces L.name with gene if L.name and gene are different, and the gene name has no comma.", {
    expect_equal(resolve_names(pbmc_dom_built_tiny, "integrin_aMb2_complex"), "integrin_aMb2_complex")
    expect_equal(resolve_names(pbmc_dom_built_tiny, "IGF1"), "IFGF1")
})

# Test resolve_complexes
test_that("resolve_complexes replaces complex names with their components.", {
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "integrin_aMb2_complex"), c("ITGB2", "ITGAM"))
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "IL7_receptor"), c("IL7R", "IL2RG"))
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "IGF1"), "IGF1")
})

# Test avg_exp_for_complexes
test_that("function returns list with average expression of component genes in complexes.", {
    exp_mat = pbmc_dom_built_tiny@counts
    complexes = pbmc_dom_built_tiny@linkages$complexes
    expect_equal(length(avg_exp_for_complexes(exp_mat, complexes), length(complexes)))
    expect_equal(length(avg_exp_for_complexes(exp_mat, complexes)[[1]]), ncol(exp_mat))
    expect_equal(avg_exp_for_complexes(exp_mat, complexes)[[1]][1], mean(exp_mat[c("ITGB2", "ITGAM"), ]))
})

# Test get_resolved_ligands
test_that("function returns list of resolved ligands.", {
    expect_equal(get_resolved_ligands(pbmc_dom_built_tiny), c("ITGB2", "ITGAM", "IGF1", "TNF", "IL7"))
})

# Test get_ligand_expression
test_that("function returns a matrix of expression", {
    all_lig_names_resolved_lists <- get_resolved_ligands(pbmc_dom_built_tiny)
    all_lig_names_resolved <- all_lig_names_resolved_lists$lig_names
    all_lig_complexes_resolved <- all_lig_names_resolved_lists$complex_names
    lig_genes <- intersect(all_lig_names_resolved, rownames(pbmc_dom_built_tiny@counts))
    result = matrix(data = c(0, 0.0333333333333333, 0, 0.783333333333333, 0, 0.45,
        0.0333333333333333, 0.00833333333333333, 0.0166666666666667,
        0, 0.591666666666667, 0.00416666666666667, 0.4, 0, 0.00833333333333333,
        0.0166666666666667, 0.0166666666666667, 0.475, 0, 0.304166666666667,
        0.0166666666666667), nrow = 7, ncol = 3, byrow = FALSE,
        dimnames = list(c("FASLG", "TNF", "TNFSF13", "PTPRC", "integrin_aVb3_complex",
            "integrin_aMb2_complex", "IL7"),
            c("CD8_T_cell", "CD14_monocyte", "B_cell")))
    expect_equal(
        get_ligand_expression(pbmc_dom_built_tiny, 
            send_clusters = levels(pbmc_dom_built_tiny@clusters),
            lig_genes, complexes = all_lig_complexes_resolved, exp_type = "counts"), 
        result)
})

# Test get_signaling_info
test_that("function returns a data frame of signaling information", {
    all_lig_names_resolved_lists <- get_resolved_ligands(pbmc_dom_built_tiny)
    all_lig_names_resolved <- all_lig_names_resolved_lists$lig_names
    all_lig_complexes_resolved <- all_lig_names_resolved_lists$complex_names
    lig_genes <- intersect(all_lig_names_resolved, rownames(pbmc_dom_built_tiny@counts))
    cl_ligands <- get_ligand_expression(dom, send_clusters, lig_genes, all_lig_complexes_resolved, exp_type)
    cl_ligands_sub <- reshape2::melt(cl_ligands)
    colnames(cl_ligands_sub) <- c("ligand", "cluster", "mean_counts")
    result = structure(list(ligand = structure(c(4L, 4L, 4L, 5L, 6L, 5L, 6L,
    5L, 6L, 7L, 7L, 7L, 4L, 4L, 4L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    levels = c("FASLG", "TNF", "TNFSF13", "PTPRC", "integrin_aVb3_complex",
    "integrin_aMb2_complex", "IL7"), class = "factor"), 
    receptor = c(rep("CD22", 3), rep("FCER2", 6), rep("IL7_receptor", 3), rep("CD22", 3),
        rep("FAS", 9)), 
    transcription_factor = c(rep("ZNF257", 9), rep("ATF4", 6), rep("RUNX1", 9)), 
    ligand_exp = c(0.783333333333333, 0.591666666666667, 0.475, 0, 0.45,
    0.00416666666666667, 0.4, 0, 0.304166666666667, 0.0333333333333333, 
    0, 0.0166666666666667, 0.783333333333333, 0.591666666666667, 0.475, 
    0, 0.0333333333333333, 0, 0.00833333333333333, 0.0166666666666667, 
    0, 0.00833333333333333, 0.0166666666666667, 0.0166666666666667), 
    rec_exp = c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
    0.504166666666667, 0.504166666666667, 0.504166666666667, 
    0.108333333333333, 0.108333333333333, 0.108333333333333, 0.0416666666666667, 
    0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 
    0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667),
    tf_auc = c(0.114904235503725, 0.114904235503725, 0.114904235503725, 
    0.114904235503725, 0.114904235503725, 0.114904235503725, 0.114904235503725, 
    0.114904235503725, 0.114904235503725, 0.114904235503725, 0.114904235503725,
    0.114904235503725, 0.128534905503555, 0.128534905503555, 0.128534905503555,
    0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 0.0477762930569053,
    0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 
    0.0477762930569053), 
    sending_cl = structure(c(1L, 2L, 3L, 1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 
    2L, 3L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), 
    levels = c("CD8_T_cell", "CD14_monocyte", "B_cell"), class = "factor"), 
    receiving_cl = c(rep("CD8_T_cell", 9), rep("CD14_monocyte", 6), rep("B_cell", 9)),
    row.names = c(NA, -24L), class = "data.frame"))
    expect_equal(
        get_signaling_info(pbmc_dom_built_tiny,
            send_clusters = levels(pbmc_dom_built_tiny@clusters), 
            lig_genes, all_lig_complexes_resolved, exp_type = "counts"),
        result)
})

# Test network_to_df
test_that("function returns data frame of ligand, receptor, tf signaling", {
    all_lig_names_resolved_lists <- get_resolved_ligands(pbmc_dom_built_tiny)
    all_lig_names_resolved <- all_lig_names_resolved_lists$lig_names
    all_lig_complexes_resolved <- all_lig_names_resolved_lists$complex_names
    lig_genes <- intersect(all_lig_names_resolved, rownames(pbmc_dom_built_tiny@counts))
    cl_ligands <- get_ligand_expression(dom, send_clusters, lig_genes, all_lig_complexes_resolved, exp_type)
    cl_ligands_sub <- reshape2::melt(cl_ligands)
    colnames(cl_ligands_sub) <- c("ligand", "cluster", "mean_counts")
    result = structure(list(ligand = structure(c(4L, 4L, 4L, 5L, 6L, 5L, 6L,
    5L, 6L, 7L, 7L, 7L, 4L, 4L, 4L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    levels = c("FASLG", "TNF", "TNFSF13", "PTPRC", "integrin_aVb3_complex",
    "integrin_aMb2_complex", "IL7"), class = "factor"), 
    receptor = c(rep("CD22", 3), rep("FCER2", 6), rep("IL7_receptor", 3), rep("CD22", 3),
        rep("FAS", 9)), 
    transcription_factor = c(rep("ZNF257", 9), rep("ATF4", 6), rep("RUNX1", 9)), 
    ligand_exp = c(0.783333333333333, 0.591666666666667, 0.475, 0, 0.45,
    0.00416666666666667, 0.4, 0, 0.304166666666667, 0.0333333333333333, 
    0, 0.0166666666666667, 0.783333333333333, 0.591666666666667, 0.475, 
    0, 0.0333333333333333, 0, 0.00833333333333333, 0.0166666666666667, 
    0, 0.00833333333333333, 0.0166666666666667, 0.0166666666666667), 
    rec_exp = c(0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
    0.504166666666667, 0.504166666666667, 0.504166666666667, 
    0.108333333333333, 0.108333333333333, 0.108333333333333, 0.0416666666666667, 
    0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 
    0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667),
    tf_auc = c(0.114904235503725, 0.114904235503725, 0.114904235503725, 
    0.114904235503725, 0.114904235503725, 0.114904235503725, 0.114904235503725, 
    0.114904235503725, 0.114904235503725, 0.114904235503725, 0.114904235503725,
    0.114904235503725, 0.128534905503555, 0.128534905503555, 0.128534905503555,
    0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 0.0477762930569053,
    0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 0.0477762930569053, 
    0.0477762930569053), 
    sending_cl = structure(c(1L, 2L, 3L, 1L, 1L, 2L, 2L, 3L, 3L, 1L, 2L, 3L, 1L, 
    2L, 3L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), 
    levels = c("CD8_T_cell", "CD14_monocyte", "B_cell"), class = "factor"), 
    receiving_cl = c(rep("CD8_T_cell", 9), rep("CD14_monocyte", 6), rep("B_cell", 9)),
    row.names = c(NA, -24L), class = "data.frame"))
    expect_equal(network_to_df(pbmc_dom_built_tiny, 
        send_clusters = levels(pbmc_dom_built_tiny@clusters), 
        rec_clusters = levels(pbmc_dom_built_tiny), 
        exp_type = "counts"), result)
})