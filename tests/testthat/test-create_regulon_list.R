test_that("create_regulon_list_scenic runs", {
    out <- create_regulon_list_scenic(regulons_tiny)
    expect_type(out, "list")
})
