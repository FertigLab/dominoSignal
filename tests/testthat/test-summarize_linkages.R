test_that("summarize_linkages runs", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    expect_s4_class(out, "linkage_summary")
})

test_that("summarise_linkages is a synonym of summarize_linkages", {
    expect_identical(summarise_linkages, summarize_linkages)

    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    out_ize <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    out_ise <- summarise_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    expect_identical(out_ize, out_ise)
})

test_that("summarize_linkages builds the expected linkage_summary structure", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)

    expect_equal(as.character(out@subject_names), c("dom1", "dom2"))
    expect_equal(out@subject_meta, meta)
    expect_named(out@subject_linkages, c("dom1", "dom2"))
    expect_named(out@subject_linkages$dom1, levels(tiny_dom1@clusters))
    expect_named(out@subject_linkages$dom1$CD8_T_cell, c("tfs", "rec", "incoming_lig", "tfs_rec", "rec_lig"))
})

test_that("summarize_linkages defaults subject_names to the first column of subject_meta", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta)

    expect_equal(as.character(out@subject_names), c("dom1", "dom2"))
})

test_that("summarize_linkages pairs tfs_rec and rec_lig correctly for a cluster with active linkages", {
    dom_ls <- list(dom1 = tiny_dom1)
    meta <- data.frame(ID = "dom1", group = "A", stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    cd8 <- out@subject_linkages$dom1$CD8_T_cell

    expect_equal(cd8$tfs_rec, c(
        "CREM <- IL7_receptor", "CREM <- CXCR3", "ZNF431 <- CXCR3",
        "FLI1 <- CXCR3", "FLI1 <- IL7_receptor"
    ))
    expect_equal(cd8$rec_lig, c("IL7_receptor <- IL7", "CXCR3 <- CCL20"))
})

test_that("summarize_linkages returns empty tfs_rec/rec_lig, not a bogus 'NA <- NA' entry, for a cluster with no active tfs or receptors", {
    # tiny_dom3's B_cell cluster has zero active transcription factors and receptors
    dom_ls <- list(dom3 = tiny_dom3)
    meta <- data.frame(ID = "dom3", group = "A", stringsAsFactors = FALSE)

    out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID)
    b_cell <- out@subject_linkages$dom3$B_cell

    expect_equal(b_cell$tfs, character(0))
    expect_equal(b_cell$tfs_rec, character(0))
    expect_equal(b_cell$rec_lig, character(0))
})

test_that("summarize_linkages errors when no provided subject names match domino_results", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("bogus1", "bogus2"), group = c("A", "B"), stringsAsFactors = FALSE)

    expect_error(
        summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID),
        "No provided subject names match"
    )
})

test_that("summarize_linkages warns and drops subject_names not present in domino_results", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)
    meta <- data.frame(ID = c("dom1", "dom2", "bogus"), group = c("A", "B", "C"), stringsAsFactors = FALSE)

    expect_warning(
        out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID),
        "not present in domino_results"
    )
    expect_equal(as.character(out@subject_names), c("dom1", "dom2"))
    expect_equal(out@subject_meta$ID, c("dom1", "dom2"))
})

test_that("summarize_linkages warns and restricts subject_meta when domino_results has unused entries", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2, dom3 = tiny_dom3)
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    expect_warning(
        out <- summarize_linkages(domino_results = dom_ls, subject_meta = meta, subject_names = meta$ID),
        "includes results only for provided subject names"
    )
    expect_equal(as.character(out@subject_names), c("dom1", "dom2"))
    expect_named(out@subject_linkages, c("dom1", "dom2"))
})

test_that("summarize_linkages errors when domino_results is not a named list", {
    meta <- data.frame(ID = c("dom1", "dom2"), group = c("A", "B"), stringsAsFactors = FALSE)

    expect_error(
        summarize_linkages(domino_results = list(tiny_dom1, tiny_dom2), subject_meta = meta, subject_names = meta$ID)
    )
    expect_error(
        summarize_linkages(domino_results = "not a list", subject_meta = meta, subject_names = meta$ID)
    )
})

test_that("summarize_linkages errors when subject_meta is not a data.frame", {
    dom_ls <- list(dom1 = tiny_dom1, dom2 = tiny_dom2)

    expect_error(
        summarize_linkages(domino_results = dom_ls, subject_meta = c("dom1", "dom2"), subject_names = c("dom1", "dom2"))
    )
})
