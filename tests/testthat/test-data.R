test_that("package data objects load", {
    expect_no_error(data(SCENIC))
    expect_no_error(data(PBMC))
    expect_no_error(data(CellPhoneDB))
    expect_no_error(data(DominoObjects))
    expect_no_error(data(LinkageSummary))
})
