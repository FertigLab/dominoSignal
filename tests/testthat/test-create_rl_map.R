test_that("create_rl_map_cellphonedb runs", {
    out <- create_rl_map_cellphonedb(
        genes = genes_tiny,
        proteins = proteins_tiny,
        interactions = interactions_tiny,
        complexes = complexes_tiny
    )

    expect_s3_class(out, "data.frame")
    expect_gt(nrow(out), 0)
})

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
