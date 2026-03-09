test_that("gene conversion functions run", {
    conv_tbl <- data.frame(
        mm.ens = c("ENSMUSG0001", "ENSMUSG0002"),
        hs.ens = c("ENSG0001", "ENSG0002"),
        mgi = c("GeneA", "GeneB"),
        hgnc = c("GENEA", "GENEB"),
        stringsAsFactors = FALSE
    )

    expect_no_error(table_convert_genes(
        genes = c("ENSG0001", "ENSG0002"),
        from = "ENSG",
        to = "HGNC",
        conversion_table = conv_tbl
    ))

})
