test_that("domino class methods run", {
    data(DominoObjects)
    dom <- DominoObjects$built_dom_tiny

    expect_s4_class(dom, "domino")
    expect_no_error(print(dom))
    expect_no_error(show(dom))
})
