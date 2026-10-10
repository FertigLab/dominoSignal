test_that("domino class methods run", {
    data(DominoObjects)
    dom <- DominoObjects$built_dom_tiny

    expect_s4_class(dom, "domino")
    expect_no_error(print(dom))
    expect_no_error(show(dom))
})

test_that("create_domino and build_domino record the package version in misc", {
    pkg_version <- as.character(packageVersion("dominoSignal"))
    dom <- create_domino(
        rl_map = rl_map_tiny,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = TRUE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    expect_equal(dom@misc$create_version, pkg_version)
    expect_null(dom@misc$build_version)

    built <- build_domino(dom, max_tf_pval = 0.05, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01)
    expect_equal(built@misc$create_version, pkg_version)
    expect_equal(built@misc$build_version, pkg_version)
})

test_that("domino validity accepts empty and existing objects", {
    expect_true(validObject(domino()))
    expect_true(validObject(tiny_created_dom1))
    expect_true(validObject(tiny_dom1))
    data(DominoObjects)
    expect_true(validObject(DominoObjects$dom_tiny))
    expect_true(validObject(DominoObjects$built_dom_tiny))
})

test_that("domino validity rejects malformed misc entries", {
    dom <- tiny_dom1
    dom@misc$build <- NA
    expect_error(validObject(dom), "misc\\$build must be a single TRUE or FALSE")

    dom <- tiny_dom1
    dom@misc$create <- c(TRUE, TRUE)
    expect_error(validObject(dom), "misc\\$create must be a single TRUE or FALSE")

    dom <- tiny_dom1
    dom@misc$create_version <- "not a version"
    expect_error(validObject(dom), "misc\\$create_version must be a single valid version string")

    dom <- tiny_dom1
    dom@misc$build_version <- c("1.0.0", "1.1.0")
    expect_error(validObject(dom), "misc\\$build_version must be a single valid version string")

    dom <- tiny_dom1
    dom@misc$build_version <- "1.7.0"
    expect_true(validObject(dom))
})

test_that("domino validity rejects matrices whose cells do not match clusters", {
    dom <- tiny_dom1
    dom@z_scores <- dom@z_scores[, -1, drop = FALSE]
    expect_error(validObject(dom), "z_scores columns must be the same cells in the same order as names\\(clusters\\)")

    dom <- tiny_dom1
    colnames(dom@counts)[1] <- "not_a_cell"
    expect_error(validObject(dom), "counts columns must be the same cells in the same order as names\\(clusters\\)")

    dom <- tiny_dom1
    dom@features <- dom@features[, -1, drop = FALSE]
    expect_error(validObject(dom), "features columns must be the same cells in the same order as names\\(clusters\\)")

    # cells are indexed by position, so the same cells in a different order are invalid
    dom <- tiny_dom1
    dom@z_scores <- dom@z_scores[, rev(colnames(dom@z_scores)), drop = FALSE]
    expect_error(validObject(dom), "z_scores columns must be the same cells in the same order")
})

test_that("domino validity rejects clust_de and cor that do not match clusters and features", {
    dom <- tiny_dom1
    colnames(dom@clust_de)[1] <- "not_a_cluster"
    expect_error(validObject(dom), "clust_de columns must match levels\\(clusters\\)")

    dom <- tiny_dom1
    dom@cor <- dom@cor[, -1, drop = FALSE]
    expect_error(validObject(dom), "cor columns must match features rows")
})

test_that("domino constructor runs validity checks", {
    clusters <- factor(c(c1 = "A", c2 = "B"))
    z <- matrix(0, nrow = 1, ncol = 2, dimnames = list("g1", c("c1", "c3")))
    expect_error(domino(clusters = clusters, z_scores = z), "z_scores columns must be the same cells")
})

test_that("print and show report recorded package versions", {
    dom <- tiny_dom1
    dom@misc$create_version <- "1.7.0"
    dom@misc$build_version <- "1.8.0"
    expect_message(print(dom), "Created with dominoSignal v1\\.7\\.0\nBuilt with dominoSignal v1\\.8\\.0")
    expect_message(show(dom), "Created with dominoSignal v1\\.7\\.0\nBuilt with dominoSignal v1\\.8\\.0")

    created <- tiny_created_dom1
    created@misc$create_version <- "1.7.0"
    expect_message(print(created), "has not been built\nCreated with dominoSignal v1\\.7\\.0")
    expect_message(show(created), "has not been built\nCreated with dominoSignal v1\\.7\\.0")
    expect_no_match(capture.output(show(created), type = "message"), "Built with")
})

test_that("print and show omit versions when none are recorded", {
    expect_no_match(capture.output(print(tiny_dom1), type = "message"), "dominoSignal v")
    expect_no_match(capture.output(show(tiny_dom1), type = "message"), "dominoSignal v")
    expect_no_match(capture.output(show(tiny_created_dom1), type = "message"), "dominoSignal v")
})

test_that("domino validity rejects empty cell-level matrices once other cell data is present", {
    dom <- tiny_dom1
    dom@counts <- methods::new("dgCMatrix")
    expect_error(validObject(dom), "counts columns must be the same cells in the same order as names\\(clusters\\)")

    dom <- tiny_dom1
    dom@features <- matrix(nrow = 0, ncol = 0)
    expect_error(validObject(dom), "features columns must be the same cells in the same order as names\\(clusters\\)")

    # matrices present without clusters
    dom <- tiny_dom1
    dom@clusters <- factor()
    expect_error(validObject(dom), "clusters must be a factor named by cell")
})

test_that("domino validity rejects unnamed clusters", {
    dom <- tiny_dom1
    names(dom@clusters) <- NULL
    expect_error(validObject(dom), "clusters must be a factor named by cell")
    expect_error(
        domino(clusters = factor(c("A", "B"))),
        "clusters must be a factor named by cell"
    )
})

test_that("domino validity rejects malformed creation variables", {
    dom <- tiny_dom1
    dom@misc$create_vars <- list(
        tf_selection_method = "clusters", tf_variance_quantile = 0.5,
        use_complexes = TRUE, rec_min_thresh = 0.025, remove_rec_dropout = FALSE
    )
    expect_true(validObject(dom))

    bad <- dom
    bad@misc$create_vars$tf_selection_method <- "bogus"
    expect_error(validObject(bad), "misc\\$create_vars\\$tf_selection_method is missing or invalid")
    bad <- dom
    bad@misc$create_vars$tf_variance_quantile <- 2
    expect_error(validObject(bad), "misc\\$create_vars\\$tf_variance_quantile is missing or invalid")
    bad <- dom
    bad@misc$create_vars$use_complexes <- NA
    expect_error(validObject(bad), "misc\\$create_vars\\$use_complexes is missing or invalid")
    bad <- dom
    bad@misc$create_vars$rec_min_thresh <- -1
    expect_error(validObject(bad), "misc\\$create_vars\\$rec_min_thresh is missing or invalid")
    bad <- dom
    bad@misc$create_vars$remove_rec_dropout <- NULL
    expect_error(validObject(bad), "misc\\$create_vars\\$remove_rec_dropout is missing or invalid")
})

test_that("domino validity rejects malformed build variables", {
    dom <- tiny_dom1
    expect_true(validObject(dom))

    bad <- dom
    bad@misc$build_vars["max_tf_per_clust"] <- -1
    expect_error(validObject(bad), "misc\\$build_vars\\[\"max_tf_per_clust\"\\] is missing or invalid")
    bad <- dom
    bad@misc$build_vars["rec_tf_cor_threshold"] <- 1.5
    expect_error(validObject(bad), "misc\\$build_vars\\[\"rec_tf_cor_threshold\"\\] is missing or invalid")
    bad <- dom
    bad@misc$build_vars <- bad@misc$build_vars[names(bad@misc$build_vars) != "min_rec_percentage"]
    expect_error(validObject(bad), "misc\\$build_vars\\[\"min_rec_percentage\"\\] is missing or invalid")
    bad <- dom
    bad@misc$build_vars <- as.character(bad@misc$build_vars)
    expect_error(validObject(bad), "misc\\$build_vars\\[\"max_tf_pval\"\\] is missing or invalid")

    # Inf is an accepted maximum
    ok <- dom
    ok@misc$build_vars["max_rec_per_tf"] <- Inf
    expect_true(validObject(ok))
})

test_that("domino validity accepts and checks the pre-rename min_tf_pval build variable", {
    # stored fixture was built before max_tf_pval was renamed
    expect_true("min_tf_pval" %in% names(tiny_dom1@misc$build_vars))
    expect_false("max_tf_pval" %in% names(tiny_dom1@misc$build_vars))
    expect_true(validObject(tiny_dom1))

    renamed <- tiny_dom1
    names(renamed@misc$build_vars)[names(renamed@misc$build_vars) == "min_tf_pval"] <- "max_tf_pval"
    expect_true(validObject(renamed))

    bad <- tiny_dom1
    bad@misc$build_vars["min_tf_pval"] <- -1
    expect_error(validObject(bad), "misc\\$build_vars\\[\"max_tf_pval\"\\] is missing or invalid")
})

test_that("print and show report creation variables when recorded", {
    dom <- tiny_dom1
    dom@misc$create_vars <- list(
        tf_selection_method = "variable", tf_variance_quantile = 0.5,
        use_complexes = FALSE, rec_min_thresh = 0.025, remove_rec_dropout = FALSE
    )
    expect_message(print(dom), "Created using 'variable' TF selection without receptor and ligand complexes")
    expect_message(show(dom), "Created using 'variable' TF selection without receptor and ligand complexes")

    created <- tiny_created_dom1
    created@misc$create_vars <- dom@misc$create_vars
    created@misc$create_vars$use_complexes <- TRUE
    expect_message(print(created), "has not been built\nCreated using 'variable' TF selection with receptor and ligand")

    expect_no_match(capture.output(print(tiny_dom1), type = "message"), "TF selection")
})
