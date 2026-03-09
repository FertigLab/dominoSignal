test_that("bool conversion function works", {
    df <- data.frame(
        list(
            c1 = c("True", "False"),
            c2 = c("False", "True"),
            c3 = c(1, 2)
        ),
        c4 = c("a", "b"),
        stringsAsFactors = FALSE
    )

    c_df <- conv_py_bools(df)

    expect_type(c_df$c1, "logical")
    expect_type(c_df$c2, "logical")
    expect_type(c_df$c3, "double")
    expect_type(c_df$c4, "character")
})

test_that("read if char tries to read a file", {
    # Handling the warning that follows the error with suppressWarnings()
    # Since error should prevent warning in normal circumstances (but not in test)
    expect_error(read_if_char(
        file.path(".", "file_that_not_exists.csv") %>% suppressWarnings()),
        "cannot open the connection")
        
    expect_error(read_if_char(c("a", "b")), "Length of obj must be one of: 1")
})

test_that("mandatory field absence in check_arg yields error, presence does not", {
    expect_error(check_arg(
        arg = data.frame(a = c(1, 2), b = c(3, 4)),
        allow_class = "data.frame",
        need_vars = "c"
    ))
    expect_silent(check_arg(
        arg = data.frame(a = c(1, 2), b = c(3, 4)),
        allow_class = "data.frame",
        need_vars = c("a", "b")
    ))
})

test_that("range checker in check_arg works", {
    expect_error(
        check_arg(1, allow_class = "numeric", allow_range = c(2, 5)),
        "All values in 1 must be between 2 and 5"
    )
    expect_silent(check_arg(3, allow_class = "numeric", allow_range = c(2, 5)))
})

test_that("check_arg works for class", {
    expect_error(
        check_arg(1, allow_class = "data.frame"),
        "Class of 1 must be one of: data.frame"
    )
    expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
            allow_class = "data.frame"
        ))
    expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
            allow_class = c("data.frame", "matrix")
        ))
})

test_that("resolve_complexes maps complex names to component genes, resolve_names maps ligand name to gene name", {
    data(PBMC)
    data(SCENIC)

    entry_for_resolve_names <- c(
        int_pair = "12oxoLeukotrieneB4_byPTGR1 & LTB4R",
        name_A = "12oxoLeukotrieneB4_byPTGR1", uniprot_A = "Q14914", gene_A = "PTGR1", type_A = "L",
        name_B = "LTB4R", uniprot_B = "Q15722", gene_B = "LTB4R", type_B = "R",
        annotation_strategy = "curated",
        source = "HMRbase;uniprot;reactome",
        database_name = "CellPhoneDB"
    )
    rl_map_tiny <- rbind(rl_map_tiny, entry_for_resolve_names)
    
    pbmc_dom_tiny <- create_domino(
        rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
        counts = PBMC$count_tiny, z_scores = PBMC$zscore_tiny,
        clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
        use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE
    )

    # Test resolve_names
    genes <- c("integrin_a6b4_complex", "TGFB3", "12oxoLeukotrieneB4_byPTGR1")
    ligand_names_resolved <- c(
        "integrin_a6b4_complex" = "integrin_a6b4_complex",
        "TGFB3" = "TGFB3",
        "12oxoLeukotrieneB4_byPTGR1" = "PTGR1"
    )
    expect_equal(resolve_names(pbmc_dom_tiny, genes), ligand_names_resolved)

    # Test resolve_complexes
    genes <- c("integrin_a6b4_complex", "IL7_receptor", "TGFBR3")
    complexes_resolved <- list(
        "integrin_a6b4_complex" = c("ITGB4", "ITGA6"),
        "IL7_receptor" = c("IL7R", "IL2RG"),
        "TGFBR3" = "TGFBR3"
    )
    expect_equal(resolve_complexes(pbmc_dom_tiny, genes), complexes_resolved)
})

test_that("avg_exp_for_complexes returns list with average expression of component genes for any complexes.", {
    exp_mat <- data.frame(
        CD8_T_cell = c(0.252464, 0.000000, 0.000000, 0.000000),
        CD14_monocyte = c(0.000000, 0.000000, 0.000000, 0.000000),
        B_cell = c(0.000000, 0.000000, 0.329416, 0.000000)
    )
    rownames(exp_mat) <- c("TGFBR3", "ITGA6", "ITGB4", "IL7R")
    genes_complexes_resolved <- list(
        "integrin_a6b4_complex" = c("ITGB4", "ITGA6"), # both genes present in data, output mean
        "IL7_receptor" = c("IL7R", "IL2RG"), # only one of the genes present in data - will be filtered
        "TGFBR3" = "TGFBR3", # not a complex, output as is
        "IL7R" = "IL7R"
    ) # gene can interact as part of a complex or alone, output as is
    complexes_avg <- list(
        "integrin_a6b4_complex" = data.frame(CD8_T_cell = 0, CD14_monocyte = 0, B_cell = 0.164708),
        "TGFBR3" = data.frame(CD8_T_cell = 0.252464, CD14_monocyte = 0, B_cell = 0, row.names = "TGFBR3"),
        "IL7R" = data.frame(CD8_T_cell = 0, CD14_monocyte = 0, B_cell = 0, row.names = "IL7R")
    )
    expect_equal(avg_exp_for_complexes(exp_mat, genes_complexes_resolved), complexes_avg)
})

test_that("mean_exp_by_cluster returns gene expression averaged over clusters.", {
    
    clust <- c("B_cell", "CD14_monocyte")
    genes <- c("CCL20", "IL7", "TGFB3")
    mean_exp <- data.frame(
        B_cell = c(0.000000, 0.3533981, 0.000000),
        CD14_monocyte = c(0.0000000, 0.0000000, 0.1461621),
        row.names = c("CCL20", "IL7", "TGFB3")
    )
    expect_equal(mean_exp_by_cluster(dom = tiny_dom1, clusts = clust, genes = genes), mean_exp, tolerance = 1e-06)
})

test_that("lc function runs", {
    expect_equal(lc(list(a = 1:2, b = 3:4), "a"), 1:2)
})

test_that("ggplot_col_gen function runs", {
    cols <- ggplot_col_gen(4)
    expect_type(cols, "character")
    expect_length(cols, 4)
})