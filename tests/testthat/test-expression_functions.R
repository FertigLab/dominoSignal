test_that("expression functions run", {
    dom <- tiny_dom1

    exp_mat <- as.data.frame(dom@z_scores[1:3, 1:3, drop = FALSE])
    clusts <- levels(dom@clusters)[1:2]
    genes <- rownames(dom@z_scores)[1:3]
    ligands <- genes[1:2]
    barcodes <- colnames(dom@counts)[1:10]

    expect_no_error(avg_exp_for_complexes(exp_mat, list(g1 = genes[1], g2 = genes[2:3])))
    expect_no_error(mean_exp_by_cluster(dom = dom, clusts = clusts, genes = genes))
    expect_no_error(mean_ligand_expression(
        x = as.matrix(dom@counts),
        ligands = ligands,
        cell_ident = clusts[1],
        cell_barcodes = barcodes,
        destination = "R_TEST"
    ))
    expect_no_error(do_norm(matrix(c(1, 2, 3, 4), nrow = 2), "row"))
})
