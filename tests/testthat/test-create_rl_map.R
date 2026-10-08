# (fake) test orthologs from tiny genes
tiny_hs_genes <- unique(genes_tiny$gene_name)
tiny_orthos <- data.frame(
    hgnc = tiny_hs_genes,
    mgi = paste0(substr(tiny_hs_genes, 1, 1), tolower(substring(tiny_hs_genes, 2))), stringsAsFactors = FALSE
)

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
    expect_identical(out_nrg1, expected)
})

test_that("create_rl_map_cellphonedb returns an empty but correctly shaped result along with a warning when no interactions are kept", {
    no_orthos <- data.frame(hgnc = "NOGENE", mgi = "Notagene")
    expect_warning(
        out <- suppressMessages(create_rl_map_cellphonedb(
            genes = genes_tiny, proteins = proteins_tiny, interactions = interactions_tiny,
            complexes = complexes_tiny, gene_conv = c("HGNC", "MGI"),
            alternate_convert = TRUE, alternate_convert_table = no_orthos
        )),
        "No receptor-ligand interactions were retained after filtering for orthologs."
    )
    expect_s3_class(out, "data.frame")
    expect_shape(out, dim = c(0, ncol(rl_map_tiny)))
    expect_identical(colnames(out), colnames(rl_map_tiny))
    expect_true(all(vapply(out, is.character, logical(1))))
})

test_that("create_rl_map_cellphonedb skips complex with component lacking an ortholog and names relevant gene", {
    complexes_extra <- complexes_tiny
    complexes_extra$uniprot_3[complexes_extra$complex_name == "IL7_receptor"] <- "FAKEUNI"
    fake_row <- genes_tiny[genes_tiny$gene_name == "IL2RG", , drop = FALSE]
    fake_row$gene_name <- "FAKE1"
    fake_row$uniprot <- "FAKEUNI"
    genes_extra <- rbind(genes_tiny, fake_row)
    msgs <- testthat::capture_messages(out <- create_rl_map_cellphonedb(
        genes = genes_extra, proteins = proteins_tiny, interactions = interactions_tiny,
        complexes = complexes_extra, gene_conv = c("HGNC", "MGI"),
        alternate_convert = TRUE, alternate_convert_table = tiny_orthos
    ))
    expect_true(any(grepl("No gene orthologs found for: FAKE1", msgs)))
    expect_false(any(grepl("IL2RG", msgs)))
    expect_true(any(grepl("Skipping interaction: IL7_receptor", msgs)))
    expect_false("IL7_receptor" %in% c(out$name_A, out$name_B))
    expect_shape(out, nrow = nrow(rl_map_tiny) - 1)
})

# create IL2RG replacement with two genes in given order
il2rg_genes <- function(gene1, gene2) {
    extra <- genes_tiny[genes_tiny$gene_name == "IL2RG", , drop = FALSE]
    rows <- extra[c(1, 1), , drop = FALSE]
    rows$gene_name <- c(gene1, gene2)
    new <- rbind(genes_tiny[genes_tiny$gene_name != "IL2RG", , drop = FALSE], rows)
    return(new)
}

test_that("create_rl_map_cellphonedb uses the gene with an ortholog for a multi-gene complex component", {
    # regardless of order, IL2RG has an ortholog and FAKE1 does not
    for (order in list(c("IL2RG", "FAKE1"), c("FAKE1", "IL2RG"))) {
        msgs <- testthat::capture_messages(out <- create_rl_map_cellphonedb(
            genes = il2rg_genes(order[1], order[2]), proteins = proteins_tiny, interactions = interactions_tiny,
            complexes = complexes_tiny, gene_conv = c("HGNC", "MGI"),
            alternate_convert = TRUE, alternate_convert_table = tiny_orthos
        ))
        expect_false(any(grepl("No gene orthologs found for", msgs, fixed = TRUE)))
        expect_equal(nrow(out), nrow(rl_map_tiny))
        expect_true("Il7r,Il2rg" %in% out$gene_A)
    }
})

test_that("create_rl_map_cellphonedb skips a complex when no gene for a component has an ortholog", {
    msgs <- testthat::capture_messages(out <- create_rl_map_cellphonedb(
        genes = il2rg_genes("FAKE1", "FAKE2"), proteins = proteins_tiny, interactions = interactions_tiny,
        complexes = complexes_tiny, gene_conv = c("HGNC", "MGI"),
        alternate_convert = TRUE, alternate_convert_table = tiny_orthos
    ))
    expect_true(any(grepl("FAKE1", msgs, fixed = TRUE)))
    expect_false(any(grepl("IL7R", msgs, fixed = TRUE)))
    expect_false("IL7_receptor" %in% c(out$name_A, out$name_B))
    expect_shape(out, nrow = nrow(rl_map_tiny) - 1)
})
