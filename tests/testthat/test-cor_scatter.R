test_that("cor_scatter runs", {
    dom <- tiny_dom1
    tf <- colnames(dom@cor)[1]
    rec <- rownames(dom@cor)[1]

    plt <- cor_scatter(dom = dom, tf = tf, rec = rec)
    expect_s3_class(plt, "gg")
})
