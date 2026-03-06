test_that("create_rl_map_cellphonedb fails on wrong input arg type.", {

    expect_error(create_rl_map_cellphonedb(
        genes = list(), proteins = proteins_tiny,
        interactions = interactions_tiny, complexes = complexes_tiny
    ), "Class of genes must be one of: character,data.frame")

    expect_error(create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = list(),
        interactions = interactions_tiny, complexes = complexes_tiny
    ), "Class of proteins must be one of: character,data.frame")

    expect_error(create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_tiny,
        interactions = list(), complexes = complexes_tiny
    ), "Class of interactions must be one of: character,data.frame")

    expect_error(create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_tiny,
        interactions = interactions_tiny, complexes = list()
    ), "Class of complexes must be one of: character,data.frame")

    expect_error(create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_tiny,
        interactions = interactions_tiny, complexes = complexes_tiny,
        database_name = list()
    ), "Class of database_name must be one of: character")

    expect_error(create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_tiny,
        interactions = interactions_tiny, complexes = complexes_tiny,
        database_name = c("length", ">1")
    ), "Length of database_name must be one of: 1")
})


test_that("create_domino fails on wrong input arg type.", {

    # bad rl map
    bad_rl_map <- "rl_map"
    expect_error(
        create_domino(bad_rl_map,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "Class of rl_map must be one of: data.frame"
    )

    bad_rl_map <- rl_map_tiny
    colnames(bad_rl_map) <- paste(colnames(bad_rl_map), "qq")
    expect_error(
        create_domino(bad_rl_map,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "Required variables gene_A, gene_B, type_A, type_B not found in rl_map"
    )

    # bad features
    bad_features <- matrix()
    expect_error(
        create_domino(rl_map_tiny,
            bad_features,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "No rownames found in features"
    )

    # seurat or counts, zscores and clusters
    expect_error(
        create_domino(
            rl_map_tiny,
            tiny_auc1
        ),
        "Class of counts must be one of: matrix,data.frame"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1
        ),
        "Class of clusters must be one of: factor"
    )

    # bad rec_min threshold
    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            rec_min_thresh = 20
        ),
        "All values in rec_min_thresh must be between 0 and 1"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            rec_min_thresh = -20
        ),
        "All values in rec_min_thresh must be between 0 and 1"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            tf_selection_method = "non-existent"
        ),
        "All values in tf_selection_method must be one of: clusters, variable, all"
    )
})
