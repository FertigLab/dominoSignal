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

# (fake) test orthologs from tiny genes
tiny_hs_genes <- unique(genes_tiny$gene_name)
tiny_orthos <- data.frame(
    hgnc = tiny_hs_genes,
    mgi = paste0(substr(tiny_hs_genes, 1, 1), tolower(substring(tiny_hs_genes, 2))), stringsAsFactors = FALSE
)

test_that("create_rl_map_cellphonedb converts genes, including complex components, with a conversion table", {
    out <- create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_tiny, interactions = interactions_tiny,
        complexes = complexes_tiny, gene_conv = c("HGNC", "MGI"),
        alternate_convert = TRUE, alternate_convert_table = tiny_orthos
    )
    expect_equal(nrow(out), nrow(rl_map_tiny))
    expect_true("Itgb4,Itga6" %in% c(out$gene_A, out$gene_B))
    expect_true(all(c(out$gene_A, out$gene_B) %in% c(tiny_orthos$mgi, "Itgb4,Itga6", "Il7r,Il2rg")))
})

test_that("create_rl_map_cellphonedb skips interactions with a complex component lacking an ortholog", {
    expect_message(
        out <- create_rl_map_cellphonedb(
            genes = genes_tiny, proteins = proteins_tiny, interactions = interactions_tiny,
            complexes = complexes_tiny, gene_conv = c("HGNC", "MGI"),
            alternate_convert = TRUE,
            alternate_convert_table = tiny_orthos[tiny_orthos$hgnc != "ITGA6", ]
        ),
        "No gene orthologs found for: ITGA6"
    )
    expect_equal(nrow(out), nrow(rl_map_tiny) - 1)
    expect_false(any(grepl("ITGA6|Itga6", c(out$gene_A, out$gene_B))))
    expect_false("integrin_a6b4_complex" %in% c(out$name_A, out$name_B))
})

test_that("create_rl_map_cellphonedb handles missing receptor annotations in proteins input", {
    proteins_na <- proteins_tiny
    proteins_na$receptor[1] <- NA
    expect_no_error(out <- create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_na, interactions = interactions_tiny, complexes = complexes_tiny
    ))
    expect_identical(out, rl_map_tiny)

    proteins_nrg1 <- proteins_tiny
    proteins_nrg1$receptor[proteins_nrg1$protein_name == "NRG1_HUMAN"] <- NA
    expect_no_error(out_nrg1 <- create_rl_map_cellphonedb(
        genes = genes_tiny, proteins = proteins_nrg1, interactions = interactions_tiny, complexes = complexes_tiny
    ))
    expected <- rl_map_tiny[!(rl_map_tiny$gene_A == "NRG1" | rl_map_tiny$gene_B == "NRG1"), ]
    identical(out_nrg1, expected)
})

