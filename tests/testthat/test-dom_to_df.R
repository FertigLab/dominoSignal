test_that("dom_to_df function runs", {
    expect_no_error(dom_to_df(tiny_dom1, exp_type = "counts"))
    expect_no_error(dom_to_df(tiny_dom2, exp_type = "z_scores"))
})

test_that("get_resolved_ligands function runs", {
    expect_no_error(get_resolved_ligands(tiny_dom1))
    expect_named(get_resolved_ligands(tiny_dom1), c("lig_names", "complex_names"))
})

test_that("get_ligand_expression function runs", {
    lig_names <- c("ITGB4", "ITGA6", "CCL20")
    complex_names <- list(integrin_a6b4_complex = c("ITGB4", "ITGA6"), CCL20 = "CCL20")
    expect_no_error(get_ligand_expression(tiny_dom1, levels(tiny_dom1@clusters), lig_names,
            complex_names, exp_type = "counts"))
})

test_that("get_signaling_info function runs", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    expect_no_error(get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
            cl_ligands_sub = ligs, exp_type = "counts"))
    expect_no_error(get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
            cl_ligands_sub = ligs, exp_type = "z_scores"))
})

test_that("get_resolved_ligands handles duplicate aliases that resolve to one gene", {
    test_dom <- dom_dup_ligand_aliases(
        tiny_dom1,
        alias_names = c("CCL20_alias_1", "CCL20_alias_2"),
        target_gene = "CCL20",
        receptor_slots = c("CCR6", "CXCR3")
    )

    expect_no_error(lig_info <- get_resolved_ligands(test_dom))
    expect_type(lig_info$lig_names, "character")
    expect_length(lig_info$lig_names, length(unique(lig_info$lig_names)))
    expect_gt(anyDuplicated(names(lig_info$complex_names)), 0)
    expect_true(all(unlist(lig_info$complex_names) %in% lig_info$lig_names))
})

test_that("get_resolved_ligands returns unique ligands and complexes", {
    test_dom <- dom_dup_ligand_aliases(
        tiny_dom1,
        alias_names = c("CCL20_alias_1", "CCL20_alias_2"),
        target_gene = "CCL20",
        receptor_slots = c("CCR6", "CXCR3")
    )

    ligs <- get_resolved_ligands(tiny_dom1)
    expect_type(ligs$lig_names, "character")
    expect_length(ligs$lig_names, length(unique(ligs$lig_names)))
    expect_true(all(lengths(ligs$complex_names) > 0))
    expect_true(all(unlist(ligs$complex_names) %in% ligs$lig_names))

    dup_ligs <- get_resolved_ligands(test_dom)
    expect_type(dup_ligs$lig_names, "character")
    expect_length(dup_ligs$lig_names, length(unique(dup_ligs$lig_names)))
    expect_true(all(lengths(dup_ligs$complex_names) > 0))
    expect_true(all(unlist(dup_ligs$complex_names) %in% dup_ligs$lig_names))
})

test_that("get_resolved_ligands handles objects with no complexes", {
    test_dom <- tiny_dom3
    test_dom@linkages$complexes <- list()
    expect_no_error(resolved_ligs <- get_resolved_ligands(test_dom))
    expect_length(resolved_ligs, 2)
    expect_named(resolved_ligs$complex_names, resolved_ligs$lig_names)
})

test_that("get_ligand_expression returns ligand-first data frame with expected clusters", {
    lig_names <- c("ITGB4", "ITGA6", "CCL20")
    complex_names <- list(integrin_a6b4_complex = c("ITGB4", "ITGA6"), CCL20 = "CCL20")
    clusters <- levels(tiny_dom1@clusters)[1:2]
    single_clust <- levels(tiny_dom1@clusters)[1]
    ligand_exp <- get_ligand_expression(tiny_dom1, clusters, lig_names, complex_names, exp_type = "counts")
    ligand_exp_single <- get_ligand_expression(tiny_dom1, single_clust, lig_names, complex_names, exp_type = "counts")

    expect_s3_class(ligand_exp, "data.frame")
    expect_s3_class(ligand_exp_single, "data.frame")
    expect_true("ligand" %in% colnames(ligand_exp))
    expect_true("ligand" %in% colnames(ligand_exp_single))
    expect_equal(setdiff(colnames(ligand_exp), "ligand"), clusters)
    expect_equal(setdiff(colnames(ligand_exp_single), "ligand"), single_clust)
    expect_true(all(ligand_exp$ligand %in% c(lig_names, names(complex_names))))
    expect_true(all(ligand_exp_single$ligand %in% c(lig_names, names(complex_names))))
})

test_that("get_ligand_expression handles cluster levels with no cells", {
    test_dom <- tiny_dom1
    test_dom@clusters <- factor(test_dom@clusters, levels = c(levels(test_dom@clusters), "empty_cluster"))
    lig_names <- c("ITGB4", "ITGA6", "CCL20")
    complex_names <- list(integrin_a6b4_complex = c("ITGB4", "ITGA6"), CCL20 = "CCL20")

    expect_error(get_ligand_expression(test_dom, "empty_cluster", lig_names, complex_names,
            exp_type = "z_scores"))

    ligand_empty_plus <- get_ligand_expression(test_dom, c("empty_cluster", "B_cell"),
        lig_names, complex_names, exp_type = "counts")
    expect_s3_class(ligand_empty_plus, "data.frame")
    expect_true("ligand" %in% colnames(ligand_empty_plus))
    expect_true("B_cell" %in% colnames(ligand_empty_plus))
    expect_false(all(is.na(ligand_empty_plus$B_cell)))

})

test_that("get_ligand_expression handles extra ligands", {
    lig_names <- c("ITGB4", "ITGA6", "CCL20", "EXTRA_LIGAND")
    complex_names <- list(integrin_a6b4_complex = c("ITGB4", "ITGA6"), CCL20 = "CCL20", EXTRA_LIGAND = "EXTRA_LIGAND")
    ligand_exp <- get_ligand_expression(tiny_dom1, levels(tiny_dom1@clusters)[1:2], lig_names, complex_names,
        exp_type = "counts")
    expect_s3_class(ligand_exp, "data.frame")
    expect_true("ligand" %in% colnames(ligand_exp))
    expect_equal(setdiff(colnames(ligand_exp), "ligand"), levels(tiny_dom1@clusters)[1:2])
    expect_true(all(ligand_exp$ligand %in% c(lig_names, names(complex_names))))
    expect_message(get_ligand_expression(tiny_dom1, levels(tiny_dom1@clusters)[1:2], lig_names, complex_names,
            exp_type = "counts"), "Some ligands not found in expression matrix: EXTRA_LIGAND")
})

test_that("get_ligand_expression handles duplicated complexes", {
    lig_genes <- c("ITGB4", "ITGA6")
    dup_complexes <- list(dup_integrin = c("ITGB4", "ITGA6"), dup_integrin = c("ITGB4", "ITGA6"))

    expect_no_error(lig_df <- get_ligand_expression(
        tiny_dom1,
        send_clusters = levels(tiny_dom1@clusters),
        lig_genes = lig_genes,
        complexes = dup_complexes,
        exp_type = "counts"
    ))

    expect_s3_class(lig_df, "data.frame")
    expect_true("ligand" %in% colnames(lig_df))
    expect_gt(anyDuplicated(lig_df$ligand), 0)
})

test_that("get_ligand_expression returns appropriate values for exp_type", {
    lig_names <- c("ITGB4", "ITGA6", "CCL20")
    complex_names <- list(integrin_a6b4_complex = c("ITGB4", "ITGA6"), CCL20 = "CCL20")
    ligand_counts <- get_ligand_expression(tiny_dom1, levels(tiny_dom1@clusters)[1:2], lig_names,
        complex_names, exp_type = "counts")
    ligand_zscores <- get_ligand_expression(tiny_dom1, levels(tiny_dom1@clusters)[1:2], lig_names,
        complex_names, exp_type = "z_scores")
    expect_true(all(ligand_counts >= 0))
    expect_false(all(ligand_zscores > 0))
})

test_that("get_signaling_info returns expected format", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    signaling_info <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "counts")
    expect_s3_class(signaling_info, "data.frame")
    expect_named(signaling_info, c("ligand", "receptor", "transcription_factor", "ligand_exp", "rec_exp",
            "tf_auc", "sending_cl", "receiving_cl"))
    expect_true(all(signaling_info$ligand %in% ligs$ligand))
    expect_true(all(signaling_info$sending_cl %in% ligs$cluster))
    expect_true(all(signaling_info$receiving_cl %in% ligs$cluster))
})

test_that("get_signaling_info computes correct values for rec_exp and tf_auc", {
    
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    signaling_info_counts <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "counts")
    test_row_counts <- signaling_info_counts[1, ]

    signaling_info_zscores <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "z_scores")
    test_row_zscores <- signaling_info_zscores[1, ]

    rec_sep <- unlist(resolve_complexes(tiny_dom1, test_row_counts$receptor))
    expect_rec_counts <- mean(dom_counts(tiny_dom1)[rec_sep,
            which(tiny_dom1@clusters == test_row_counts$receiving_cl)])
    expect_rec_zscores <- mean(dom_zscores(tiny_dom1)[rec_sep,
            which(tiny_dom1@clusters == test_row_zscores$receiving_cl)])
    expect_tf_auc <- mean(dom_tf_activation(tiny_dom1)[test_row_counts$transcription_factor,
            which(tiny_dom1@clusters == test_row_counts$receiving_cl)])

    expect_equal(test_row_counts$rec_exp, expect_rec_counts)
    expect_equal(test_row_zscores$rec_exp, expect_rec_zscores)
    expect_equal(test_row_counts$tf_auc, expect_tf_auc)
    expect_equal(test_row_zscores$tf_auc, expect_tf_auc)
})

test_that("get_signaling_info returns empty data frame when no interactions are found", {
    test_dom <- tiny_dom1
    test_dom@linkages$clust_tf_rec[["B_cell"]] <- list()
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    signaling_info <- get_signaling_info(test_dom, rec_clusters = "B_cell", cl_ligands_sub = ligs, exp_type = "counts")
    expect_s3_class(signaling_info, "data.frame")
    expect_length(signaling_info, 0)
    expect_message(get_signaling_info(test_dom, rec_clusters = "B_cell", cl_ligands_sub = ligs, exp_type = "counts"),
        "No interactions found for the specified clusters and expression type.")
})

test_that("get_signaling_info handles single and plural parameters", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    expect_no_error(get_signaling_info(tiny_dom1, c("B_cell", "CD8_T_cell"), cl_ligands_sub = ligs,
            exp_type = "counts"))
    expect_no_error(get_signaling_info(tiny_dom1, c("B_cell", "CD14_monocyte"), cl_ligands_sub = ligs[1, ],
            exp_type = "z_scores"))
})

test_that("get_signaling_info only returns signaling for provided ligands", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    ligs <- ligs[ligs$ligand == "CCL20", ]
    signaling_info <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "counts")
    expect_true(all(signaling_info$ligand == "CCL20"))
    expect_true(all(signaling_info$sending_cl %in% ligs$cluster))
    expect_true(all(signaling_info$receiving_cl %in% c("CD8_T_cell", "B_cell")))
})

test_that("get_signaling_info handles ligands in cl_ligands_sub that are not in the object", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "EXTRA_LIGAND"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    signaling_info <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "counts")
    expect_true(all(signaling_info$ligand == "integrin_a6b4_complex"))
})

test_that("get_signaling_info returns appropriate values for exp_type", {
    ligs <- data.frame(ligand = rep(c("integrin_a6b4_complex", "CCL20"), n = 3),
        cluster = rep(c("CD8_T_cell", "CD14_monocyte", "B_cell"), each = 2),
        mean_counts = c(0, 0.01, 0, 0, 0.0045, 0), stringsAsFactors = FALSE)
    signaling_counts <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "counts")
    signaling_zscores <- get_signaling_info(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        cl_ligands_sub = ligs, exp_type = "z_scores")
    conserved_cols <- c("ligand", "receptor", "transcription_factor", "sending_cl", "receiving_cl", "ligand_exp")
    expect_false(any(signaling_counts$rec_exp <= 0))
    expect_false(all(signaling_zscores$rec_exp > 0))
    expect_equal(signaling_counts[ , conserved_cols], signaling_zscores[ , conserved_cols])
})

test_that("dom_to_df handles duplicate resolved ligand names without rowname errors", {
    test_dom <- dom_dup_ligand_aliases(
        tiny_dom1,
        alias_names = c("CCL20_alias_1", "CCL20_alias_2"),
        target_gene = "CCL20",
        receptor_slots = c("CCR6", "CXCR3")
    )

    expect_no_error(df <- dom_to_df(test_dom, exp_type = "counts"))
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0)
})

test_that("dom_to_df handles duplicate resolved receptor names without rowname errors", {
    test_dom <- dom_dup_receptor_aliases(
        tiny_dom1,
        alias_names = c("CCR6_alias_1", "CCR6_alias_2"),
        target_gene = "CCR6",
        ligand_slots = c("CCL20", "CCL20")
    )

    expect_no_error(df <- dom_to_df(test_dom, exp_type = "counts"))
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0)
})

test_that("dom_to_df filters by send_clusters and rec_clusters and handles single or multiple values", {
    dom_df_to_1 <- dom_to_df(tiny_dom1, rec_clusters = tiny_dom1@clusters[1],
        send_clusters = levels(tiny_dom1@clusters), exp_type = "counts")
    dom_df_from_1 <- dom_to_df(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters),
        send_clusters = tiny_dom1@clusters[1], exp_type = "counts")
    dom_df_1_to_1 <- dom_to_df(tiny_dom1, rec_clusters = tiny_dom1@clusters[1],
        send_clusters = tiny_dom1@clusters[1], exp_type = "counts")
    dom_df_2_clusts <- dom_to_df(tiny_dom1, rec_clusters = levels(tiny_dom1@clusters)[1:2],
        send_clusters = levels(tiny_dom1@clusters)[2:3], exp_type = "counts")
    
    expect_true(all(dom_df_to_1$receiving_cl == levels(tiny_dom1@clusters)[1]))
    expect_true(all(dom_df_from_1$sending_cl == levels(tiny_dom1@clusters)[1]))
    expect_true(all(dom_df_1_to_1$receiving_cl == levels(tiny_dom1@clusters)[1]))
    expect_true(all(dom_df_1_to_1$sending_cl == levels(tiny_dom1@clusters)[1]))
    expect_true(all(dom_df_2_clusts$receiving_cl %in% c(levels(tiny_dom1@clusters)[1], levels(tiny_dom1@clusters)[2])))
    expect_true(all(dom_df_2_clusts$sending_cl %in% c(levels(tiny_dom1@clusters)[2], levels(tiny_dom1@clusters)[3])))
})

test_that("dom_to_df handles duplicate cluster names with informative message", {
    expect_message(dom_to_df(tiny_dom1, rec_clusters = c("B_cell", "B_cell"), exp_type = "counts"),
        "Duplicate entries in rec_clusters detected. Using unique values: B_cell")
    expect_message(dom_to_df(tiny_dom1, send_clusters = c("B_cell", "B_cell", "CD8_T_cell"), exp_type = "counts"),
        "Duplicate entries in send_clusters detected. Using unique values: B_cell, CD8_T_cell")
    df_dup_send <- dom_to_df(tiny_dom1, send_clusters = c("B_cell", "B_cell", "CD8_T_cell"),
        rec_clusters = "CD14_monocyte", exp_type = "counts")
    df_dup_rec <- dom_to_df(tiny_dom1, send_clusters = "B_cell",
        rec_clusters = c("CD14_monocyte", "CD14_monocyte", "CD8_T_cell"), exp_type = "counts")
    expect_s3_class(df_dup_send, "data.frame")
    expect_s3_class(df_dup_rec, "data.frame")
    expect_true(all(df_dup_send$sending_cl %in% c("B_cell", "CD8_T_cell")))
    expect_true(all(df_dup_send$receiving_cl == "CD14_monocyte"))
    expect_true(all(df_dup_rec$sending_cl == "B_cell"))
    expect_true(all(df_dup_rec$receiving_cl %in% c("CD14_monocyte", "CD8_T_cell")))
})

test_that("dom_to_df doesn't return rows with no signaling for expression type 'counts'", {
    dom_df <- dom_to_df(tiny_dom1, exp_type = "counts")
    dom_df_z <- dom_to_df(tiny_dom1, exp_type = "z_scores")
    dom_df_z_gtz <- dplyr::filter(dom_df_z, ligand_exp > 0, rec_exp > 0, tf_auc > 0)
    expect_true(all(dom_df$ligand_exp > 0 & dom_df$rec_exp > 0 & dom_df$tf_auc > 0))
    expect_gte(nrow(dom_df_z), nrow(dom_df_z_gtz))

    no_fli1_dom <- dom_add_zero_auc(tiny_dom1, "FLI1")
    no_fli_df_counts <- dom_to_df(no_fli1_dom, exp_type = "counts")
    no_fli_df_z <- dom_to_df(no_fli1_dom, exp_type = "z_scores")
    expect_true(all(no_fli_df_counts$ligand_exp > 0 & no_fli_df_counts$rec_exp > 0 & no_fli_df_counts$tf_auc > 0))
    expect_false(any(no_fli_df_z$tf_auc <= 0))

    no_ccr6_dom <- dom_add_zero_exp(tiny_dom1, "CCR6")
    no_ccr6_df_counts <- dom_to_df(no_ccr6_dom, exp_type = "counts")
    no_ccr6_df_z <- dom_to_df(no_ccr6_dom, exp_type = "z_scores")
    expect_true(all(no_ccr6_df_counts$ligand_exp > 0 & no_ccr6_df_counts$rec_exp > 0 & no_ccr6_df_counts$tf_auc > 0))
    expect_true(all(dplyr::filter(no_ccr6_df_z, receptor == "CCR6")$rec_exp == 0))
})

test_that("dom_to_df returns all clusters when send_clusters and rec_clusters are not provided", {
    dom_df_all <- dom_to_df(tiny_dom1, exp_type = "counts")
    expect_gt(sum(levels(tiny_dom1@clusters) %in% dom_df_all$sending_cl), 1)
    expect_gt(sum(levels(tiny_dom1@clusters) %in% dom_df_all$receiving_cl), 1)
})

test_that("dom_to_df errors are informative", {
    expect_error(dom_to_df(tiny_dom1, rec_clusters = "nonexistent_cluster", exp_type = "counts"),
        "rec_clusters must be: CD8_T_cell, CD14_monocyte, B_cell")
    expect_error(dom_to_df(tiny_dom1, send_clusters = "nonexistent_cluster", exp_type = "counts"),
        "send_clusters must be: CD8_T_cell, CD14_monocyte, B_cell")
    expect_error(dom_to_df(tiny_dom1, exp_type = "invalid_exp_type"),
        "All values in exp_type must be: counts, z_scores")
    expect_error(dom_to_df(tiny_dom1, exp_type = c("counts", "z_scores")),
        "Length of exp_type must be one of: 1")
})
