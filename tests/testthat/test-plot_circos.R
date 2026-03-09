test_that("plot_circos functions run", {
    receptor <- "CXCR3"
    ligands <- "CCL20"

    expect_no_error(circos_ligand_receptor(tiny_dom1, receptor = receptor))
    expect_no_error(obtain_circos_expression(tiny_dom1, receptor = receptor, ligands = ligands))

    circos_df <- data.frame(
        origin = c("A_L1", "B_L1"),
        destination = "R1",
        mean.expression = c(0.5, 0.2),
        sender = c("A", "B"),
        ligand = c("L1", "L1"),
        receptor = "R1",
        scaled.mean.expression = c(1, 0.4),
        ligand.arc = c(0.5, 0.5),
        receptor.arc = c(1, 1),
        stringsAsFactors = FALSE
    )
    expect_no_error(render_circos_ligand_receptor(circos_df, receptor = "R1"))
})

test_that("domino object interpretation by obtain_circos_expression", {
    # start temporary graphics device for testing to preserve package enviroment
    png(filename = file.path(tempdir(), "ts.png"))

    # parent function runs to completion without error
    expect_no_error(circos_ligand_receptor(tiny_dom1, receptor = "CXCR3"))

    # specification of a single ligand
    circos_df_ccl20 <- obtain_circos_expression(tiny_dom1, receptor = "CXCR3", ligands = "CCL20")
    expect_equal(unique(circos_df_ccl20$ligand), "CCL20")

    # fail without specification of a ligands
    expect_error(obtain_circos_expression(tiny_dom1, receptor = "CXCR3"), "argument \"ligands\" is missing, with no default")

    dev.off()
})

test_that("Cell types with hyphenated names can be plotted", {
    # start temporary graphics device for testing to preserve package enviroment
    png(filename = file.path(tempdir(), "ts.png"))

    ts_signaling_df <- data.frame(
        origin = c("CT_1-L1", "CT_1-L2", "CT_2-L1", "CT_2-L2", "CT_3-L1", "CT_3-L2"),
        destination = "R1",
        mean.expression = c(0, 0.5, 1, 0.5, 0.5, 0.5),
        sender = c("CT_1", "CT_1", "CT_2", "CT_2", "CT_3", "CT_3"),
        ligand = c("L1", "L2", "L1", "L2", "L1", "L2"),
        receptor = "R1",
        scaled.mean.expression = c(0, 0.5, 1, 0.5, 0.5, 0.5),
        ligand.arc = 1,
        receptor.arc = 4 / 6,
        stringsAsFactors = FALSE
    )

    # plot can render from a manually written signaling data.frame
    expect_no_error(render_circos_ligand_receptor(ts_signaling_df, receptor = "R1"))

    # plot can render if all cell type underscores are replaced by hyphens
    ts_signaling_df_hyphen <- ts_signaling_df
    ts_signaling_df_hyphen$origin <- gsub("_", "-", ts_signaling_df$origin, fixed = TRUE)
    ts_signaling_df_hyphen$sender <- gsub("_", "-", ts_signaling_df$sender, fixed = TRUE)

    expect_no_error(render_circos_ligand_receptor(ts_signaling_df_hyphen, receptor = "R1"))

    dev.off()
})

test_that("Plots for receptors that have ligands in the rl_map but not the signaling matrix do not fail", {
    # Manually add ligand receptor interactions not presenting signaling
    rl_map_append <- rbind(
        rl_map_tiny,
        data.frame(
            "int_pair" = c("CXCR3 & CCX", "CXCR3 & CCY"),
            "name_A" = c("CXCR3", "CXCR3"),
            "uniprot_A" = c("P49682", "P49682"),
            "gene_A" = c("CXCR3", "CXCR3"),
            "type_A" = c("R", "R"),
            "name_B" = c("CCX", "CCY"),
            "uniprot_B" = c("CCX", "CCY"),
            "gene_B" = c("CCX", "CCY"),
            "type_B" = c("L", "L"),
            "annotation_strategy" = c("unit_test", "unit_test"),
            "source" = c("unit_test", "unit_test"),
            "database_name" = c("unit_test", "unit_test"),
            stringsAsFactors = FALSE
        )
    )

    pbmc_dom_tiny <- create_domino(
        rl_map = rl_map_append, features = tiny_auc2,
        counts = tiny_counts2, z_scores = tiny_zscores2,
        clusters = tiny_clusters2, tf_targets = regulon_list_tiny,
        use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE
    )

    dom <- build_domino(
        dom = pbmc_dom_tiny, min_tf_pval = 0.05, max_tf_per_clust = Inf,
        max_rec_per_tf = Inf, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01
    )
    # Need to adapt to testthat(3) handling of messages at some point
    testthat::local_edition(2)

    expect_message(
        obtain_circos_expression(dom = dom, receptor = "CXCR3", ligands = c("CCL20", "CCX", "CCY")),
        "Ligands: CCX, CCY of receptor CXCR3 are listed in the rl_map, but not present in the signaling matrix.",
        "Only ligands: CCL20 will be considered."
    )

    # The same message is returned from the full function creating the circos plot
    expect_message(
        circos_ligand_receptor(dom = dom, receptor = "CXCR3"),
        "Ligands: CCX, CCY of receptor CXCR3 are listed in the rl_map, but not present in the signaling matrix.",
        "Only ligands: CCL20 will be considered."
    )

    # A receptor with no ligands present returns an error
    expect_error(
        obtain_circos_expression(dom = dom, receptor = "CXCR3", ligands = c("CCX", "CCY")),
        "No ligands of receptor CXCR3 present in signaling matrix."
    )

})
