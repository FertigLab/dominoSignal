test_that("summarize_linkages runs", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    expect_s4_class(out, "linkage_summary")
})
