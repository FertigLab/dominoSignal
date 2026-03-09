test_that("linkage_summary class methods run", {
    expect_s4_class(tiny_linkage_summary, "linkage_summary")
    expect_no_error(print(tiny_linkage_summary))
    expect_no_error(show(tiny_linkage_summary))
})
