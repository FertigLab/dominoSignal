test_that("create_domino runs with tiny inputs", {
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

    expect_s4_class(dom, "domino")
    # stored fixture predates version and creation variable recording, so compare without them
    expect_equal(dom@misc$create_version, as.character(packageVersion("dominoSignal")))
    expect_equal(dom@misc$create_vars, list(
        tf_selection_method = "clusters", tf_variance_quantile = 0.5,
        use_complexes = TRUE, rec_min_thresh = 0.025, remove_rec_dropout = FALSE
    ))
    dom@misc$create_version <- NULL
    dom@misc$create_vars <- NULL
    expect_equal(dom, tiny_created_dom1)
})

test_that("create_domino fails on wrong input arg type", {

    # bad rl map
    bad_rl_map <- "rl_map"
    expect_error(
        create_domino(bad_rl_map,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "Class of rl_map must be one of: data.frame"
    )

    bad_rl_map <- rl_map_tiny
    colnames(bad_rl_map) <- paste(colnames(bad_rl_map), "qq")
    expect_error(
        create_domino(bad_rl_map,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "Required variables gene_A, gene_B, type_A, type_B not found in rl_map"
    )

    # bad features
    bad_features <- matrix()
    expect_error(
        create_domino(rl_map_tiny,
            bad_features,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1
        ),
        "No rownames found in features"
    )

    # seurat or counts, zscores and clusters
    expect_error(
        create_domino(
            rl_map_tiny,
            tiny_auc1
        ),
        "missing, with no default"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1
        ),
        "missing, with no default"
    )

    # bad rec_min threshold
    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            rec_min_thresh = 20
        ),
        "All values in rec_min_thresh must be between 0 and 1"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            rec_min_thresh = -20
        ),
        "All values in rec_min_thresh must be between 0 and 1"
    )

    expect_error(
        create_domino(rl_map_tiny,
            tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            tf_selection_method = "non-existent"
        ),
        "All values in tf_selection_method must be: clusters, variable, all"
    )
})

test_that("create_domino runs using custom rl_map with only gene and type columns", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]
    dom_custom <- create_domino(
        rl_map = rl_map_custom,
        features = tiny_auc1,
        counts = tiny_counts1,
        z_scores = tiny_zscores1,
        clusters = tiny_clusters1,
        tf_targets = regulon_list_tiny,
        use_complexes = TRUE,
        remove_rec_dropout = FALSE,
        verbose = FALSE
    )
    expect_length(dom_linkages(dom_custom, "receptor-ligand"), nrow(rl_map_custom))
    expect_identical(unname(unlist(dom_linkages(dom_custom))), unname(unlist(dom_linkages(tiny_dom1))))
})

# create and build domino pipeline with and without name columns in rl_map should give
# identical signaling results when no complexes are used
test_that("create_domino and build_domino give identical signaling with and without name columns", {
    rl_map_custom <- rl_map_tiny[, c("gene_A", "gene_B", "type_A", "type_B")]

    test_doms <- function(rl_map) {
        dom <- create_domino(
            rl_map = rl_map,
            features = tiny_auc1,
            counts = tiny_counts1,
            z_scores = tiny_zscores1,
            clusters = tiny_clusters1,
            tf_targets = regulon_list_tiny,
            use_complexes = FALSE,
            remove_rec_dropout = FALSE,
            verbose = FALSE
        )

        dom <- build_domino(dom,
            max_tf_pval = 0.05,
            max_tf_per_clust = Inf,
            max_rec_per_tf = Inf,
            rec_tf_cor_threshold = 0.1,
            min_rec_percentage = 0.01
        )
        return(dom)
    }
    dom <- test_doms(rl_map_tiny)
    dom_custom <- test_doms(rl_map_custom)

    expect_identical(unname(as.matrix(dom_signaling(dom))), unname(as.matrix(dom_signaling(dom_custom))))
})

tiny_create_args <- list(
    rl_map = rl_map_tiny, features = tiny_auc1, counts = tiny_counts1, z_scores = tiny_zscores1,
    clusters = tiny_clusters1, tf_targets = regulon_list_tiny, use_complexes = TRUE,
    remove_rec_dropout = FALSE, verbose = FALSE
)

test_that("create_domino aligns clusters, counts, and features to z_scores cell order", {
    set.seed(1)
    shuffled_clusters <- tiny_clusters1[sample(length(tiny_clusters1))]
    shuffled_counts <- tiny_counts1[, sample(ncol(tiny_counts1)), drop = FALSE]
    shuffled_features <- tiny_auc1[, sample(ncol(tiny_auc1)), drop = FALSE]
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$clusters <- shuffled_clusters
    tiny_mod_args$counts <- shuffled_counts
    tiny_mod_args$features <- shuffled_features
    dom <- expect_no_warning(do.call(create_domino, tiny_mod_args))
    expect_named(dom@clusters, colnames(dom@z_scores))
    expect_identical(colnames(dom@counts), colnames(dom@z_scores))
    expect_identical(colnames(dom@features), colnames(dom@z_scores))
    expect_equal(dom, do.call(create_domino, tiny_create_args))
})

test_that("create_domino subsets to shared cells with a warning", {
    cells <- colnames(tiny_zscores1)
    keep <- cells[-(1:3)]
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$clusters <- tiny_clusters1[cells[-1]]
    tiny_mod_args$counts <- tiny_counts1[, cells[-2], drop = FALSE]
    tiny_mod_args$features <- tiny_auc1[, cells[-3], drop = FALSE]
    expect_warning(
        dom <- do.call(create_domino, tiny_mod_args),
        "Using the 327 cells shared .* z_scores = 3, counts = 2, features = 2, clusters = 2"
    )
    expect_named(dom@clusters, keep)
    expect_identical(colnames(dom@counts), keep)
    expect_identical(colnames(dom@z_scores), keep)
    expect_identical(colnames(dom@features), keep)
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$counts <- tiny_counts1[ , keep, drop = FALSE]
    tiny_mod_args$z_scores <- tiny_zscores1[ , keep, drop = FALSE]
    tiny_mod_args$features <- tiny_auc1[ , keep, drop = FALSE]
    tiny_mod_args$clusters <- tiny_clusters1[keep]
    expect_equal(dom, do.call(create_domino, tiny_mod_args))
})

test_that("create_domino errors when no cells are shared", {
    no_shared <- tiny_clusters1
    names(no_shared) <- paste0("other_", names(no_shared))
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$clusters <- no_shared
    expect_error(
        do.call(create_domino, tiny_mod_args),
        "No cells are shared"
    )
})

test_that("create_domino drops levels of clusters whose cells are all removed", {
    dropped_clust <- levels(tiny_clusters1)[1]
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$clusters <- tiny_clusters1[tiny_clusters1 != dropped_clust]
    expect_warning(
        dom <- do.call(create_domino, tiny_mod_args), # nolint: implicit_assignment_linter.
        "Using the .* cells shared"
    )
    expect_identical(levels(dom@clusters), levels(tiny_clusters1)[-1])
    expect_identical(colnames(dom@clust_de), levels(tiny_clusters1)[-1])
    expect_identical(colnames(dom@misc$cl_rec_percent), levels(tiny_clusters1)[-1])
    built <- build_domino(dom, max_tf_pval = 0.05, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01)
    expect_named(built@linkages$clust_tf, levels(tiny_clusters1)[-1])
})

test_that("variable and all TF selection methods rank TFs by correlation with expressed receptors", {
    for (method in c("variable", "all")) {
        tiny_mod_args <- tiny_create_args
        tiny_mod_args$tf_selection_method <- method
        dom <- do.call(create_domino, tiny_mod_args)
        # no differential testing, but receptor expression by cluster is still calculated
        expect_identical(dom@clusters, tiny_clusters1)
        expect_equal(dom@misc$create_vars$tf_selection_method, method)
        expect_identical(dim(dom@clust_de), c(0L, 0L))
        expect_identical(colnames(dom@misc$cl_rec_percent), levels(tiny_clusters1))

        # with no cap on receptors per TF, clust_tf_rec holds every expressed receptor above the threshold,
        # so each active TF must have one and TFs must be ordered by their maximum correlation among them
        built <- build_domino(dom,
            max_tf_per_clust = Inf, max_rec_per_tf = Inf, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01
        )
        expect_named(built@linkages$clust_tf, levels(tiny_clusters1))
        for (clust in levels(tiny_clusters1)) {
            tfs <- built@linkages$clust_tf[[clust]]
            expect_gt(length(tfs), 0)
            means <- rowMeans(dom@features[tfs, tiny_clusters1 == clust, drop = FALSE])
            expect_true(all(means > 0))
            max_cor <- vapply(tfs, function(tf) {
                recs <- built@linkages$clust_tf_rec[[clust]][[tf]]
                if (length(recs) == 0) NA_real_ else max(dom@cor[recs, tf])
            }, numeric(1))
            expect_false(anyNA(max_cor))
            expect_true(all(max_cor > 0.1))
            expect_false(is.unsorted(rev(max_cor)))
        }
        expect_named(built@linkages$clust_rec, levels(tiny_clusters1))
        expect_identical(dim(built@signaling), rep(nlevels(tiny_clusters1), 2))
        expect_named(built@cl_signaling_matrices, levels(tiny_clusters1))

        # max_tf_per_clust keeps the top ranked TFs
        capped <- build_domino(dom,
            max_tf_per_clust = 2, max_rec_per_tf = Inf, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01
        )
        expect_identical(capped@linkages$clust_tf, lapply(built@linkages$clust_tf, head, 2))
    }
})

test_that("variable and all TF selection methods exclude TFs without expressed receptors above the threshold", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$tf_selection_method <- "all"
    dom <- do.call(create_domino, tiny_mod_args)
    no_cor <- build_domino(dom, max_tf_per_clust = Inf, rec_tf_cor_threshold = 0.99, min_rec_percentage = 0.01)
    expect_true(all(lengths(no_cor@linkages$clust_tf) == 0))
    no_expr <- build_domino(dom, max_tf_per_clust = Inf, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.99)
    expect_true(all(lengths(no_expr@linkages$clust_tf) == 0))
})

test_that("variable and all TF selection methods ignore max_tf_pval and exclude TFs with no positive mean score", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$tf_selection_method <- "all"
    dom <- do.call(create_domino, tiny_mod_args)
    strict <- build_domino(dom, max_tf_pval = 0, max_tf_per_clust = Inf, rec_tf_cor_threshold = 0.1)
    loose <- build_domino(dom, max_tf_pval = 1, max_tf_per_clust = Inf, rec_tf_cor_threshold = 0.1)
    expect_identical(strict@linkages$clust_tf, loose@linkages$clust_tf)

    # a TF selected in a cluster is excluded once its scores there are zero
    clust <- levels(tiny_clusters1)[1]
    inactive_tf <- loose@linkages$clust_tf[[clust]][1]
    expect_false(is.na(inactive_tf))
    dom@features[inactive_tf, tiny_clusters1 == clust] <- 0
    built <- build_domino(dom, max_tf_per_clust = Inf, rec_tf_cor_threshold = 0.1)
    expect_false(inactive_tf %in% built@linkages$clust_tf[[clust]])
})

test_that("build_domino infers the TF selection method for objects without creation variables", {
    # cluster route when clust_de is present
    expect_null(tiny_created_dom1@misc$create_vars)
    built <- build_domino(tiny_created_dom1, max_tf_pval = 0.05, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.01)
    expect_named(built@linkages$clust_tf, levels(tiny_clusters1))

    # mean score route when clust_de is empty
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$tf_selection_method <- "all"
    dom <- do.call(create_domino, tiny_mod_args)
    expected <- build_domino(dom, rec_tf_cor_threshold = 0.1)
    dom@misc$create_vars <- NULL
    built <- build_domino(dom, rec_tf_cor_threshold = 0.1)
    expect_identical(built@linkages$clust_tf, expected@linkages$clust_tf)
})

test_that("create_domino converts dense counts to dgCMatrix with identical results", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$counts <- as.matrix(tiny_counts1)
    dom_dense <- do.call(create_domino, tiny_mod_args)
    expect_s4_class(dom_dense@counts, "dgCMatrix")
    expect_equal(dom_dense, do.call(create_domino, tiny_create_args))

    tiny_mod_args$counts <- as.data.frame(as.matrix(tiny_counts1))
    expect_equal(do.call(create_domino, tiny_mod_args), do.call(create_domino, tiny_create_args))
})

test_that("create_domino excludes receptors missing from z_scores with a warning", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$z_scores <- tiny_zscores1[rownames(tiny_zscores1) != "CXCR3", , drop = FALSE]
    expect_warning(
        dom <- do.call(create_domino, tiny_mod_args),
        "Receptor genes not found in z_scores are excluded from correlation calculations: CXCR3"
    )
    expect_false("CXCR3" %in% rownames(dom@misc$rec_cor))
    expect_false("CXCR3" %in% rownames(dom@misc$cl_rec_percent))
    expect_equal(unname(dom@cor["CXCR3", ]), rep(0, nrow(dom@features)))
})

test_that("create_domino keeps receptor names that are not syntactic R names", {
    rename_gene <- function(x) {
        rownames(x)[rownames(x) == "CXCR3"] <- "CXCR-3"
        return(x)
    }
    rl_map_dash <- rl_map_tiny
    rl_map_dash[rl_map_dash == "CXCR3"] <- "CXCR-3"
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$rl_map <- rl_map_dash
    tiny_mod_args$counts <- rename_gene(tiny_counts1)
    tiny_mod_args$z_scores <- rename_gene(tiny_zscores1)
    dom <- do.call(create_domino, tiny_mod_args)
    expect_identical(rownames(dom@cor), names(dom@linkages$rec_lig))
    expect_equal(unname(dom@cor["CXCR-3", ]), unname(tiny_created_dom1@cor["CXCR3", ]))

    # the renamed receptor must be linked in the same clusters as the original
    build_args <- list(
        max_tf_pval = 0.05, max_tf_per_clust = Inf, max_rec_per_tf = Inf, rec_tf_cor_threshold = 0.1,
        min_rec_percentage = 0.01
    )
    built <- do.call(build_domino, c(list(dom = dom), build_args))
    built_ref <- do.call(build_domino, c(list(dom = tiny_created_dom1), build_args))
    expect_identical(
        lapply(built@linkages$clust_rec, sort),
        lapply(built_ref@linkages$clust_rec, function(x) sort(sub("^CXCR3$", "CXCR-3", x)))
    )
})

test_that("create_domino keeps features above tf_variance_quantile for the variable method", {
    cv <- apply(tiny_auc1, 1, function(x) sd(x) / mean(x))
    for (q in c(0.25, 0.5, 0.75)) {
        tiny_mod_args <- tiny_create_args
        tiny_mod_args$tf_selection_method <- "variable"
        tiny_mod_args$tf_variance_quantile <- q
        dom <- do.call(create_domino, tiny_mod_args)
        expect_setequal(rownames(dom@features), names(cv)[rank(cv) > length(cv) * q])
        expect_setequal(colnames(dom@cor), rownames(dom@features))
    }

    # a single retained feature stays a one-row matrix
    tiny_mod_args$tf_variance_quantile <- 0.99
    dom <- do.call(create_domino, tiny_mod_args)
    expect_identical(dim(dom@features), c(1L, ncol(tiny_auc1)))

    tiny_mod_args$tf_variance_quantile <- 1.5
    expect_error(
        do.call(create_domino, tiny_mod_args),
        "All values in tf_variance_quantile must be between 0 and 1"
    )
})

test_that("create_domino correlates receptors with TFs only in cells with receptor counts when remove_rec_dropout = TRUE", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$remove_rec_dropout <- TRUE
    dom <- do.call(create_domino, tiny_mod_args)
    expect_true(dom@misc$create_vars$remove_rec_dropout)
    keep <- which(tiny_counts1["CXCR3", ] > 0)
    expected <- cor(tiny_zscores1["CXCR3", keep], tiny_auc1["FLI1", keep], method = "spearman")
    expect_equal(dom@misc$rec_cor["CXCR3", "FLI1"], expected)
    expect_false(isTRUE(all.equal(dom@misc$rec_cor, tiny_created_dom1@misc$rec_cor)))
})

test_that("create_domino sets correlation to NAs for receptors in a TF's regulon", {
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$tf_targets <- list(FLI1 = "CXCR3")
    dom <- do.call(create_domino, tiny_mod_args)
    expect_true(is.na(dom@misc$rec_cor["CXCR3", "FLI1"]))
    expect_gt(tiny_created_dom1@misc$rec_cor["CXCR3", "FLI1"], 0)
    # other receptor-TF pairs are unchanged
    other <- dom@misc$rec_cor
    other["CXCR3", "FLI1"] <- tiny_created_dom1@misc$rec_cor["CXCR3", "FLI1"]
    expect_equal(other, tiny_created_dom1@misc$rec_cor)
})

test_that("create_domino summarizes receptor complex correlations by the median of their components", {
    rec_cor <- tiny_created_dom1@misc$rec_cor
    expect_equal(
        tiny_created_dom1@cor["IL7_receptor", ],
        apply(rec_cor[c("IL7R", "IL2RG"), , drop = FALSE], 2, median)
    )
    expect_equal(tiny_created_dom1@cor["CXCR3", ], rec_cor["CXCR3", ])
})

test_that("create_domino excludes rl_map rows that do not pair one receptor with one ligand", {
    rl_map_bad <- rl_map_tiny
    extra <- rl_map_tiny[c(2, 3, 2), ]
    extra$type_A <- c("L", "R", NA)
    extra$type_B <- c("L", "R", "L")
    rl_map_bad <- rbind(rl_map_bad, extra)
    tiny_mod_args <- tiny_create_args
    tiny_mod_args$rl_map <- rl_map_bad
    expect_warning(
        dom <- do.call(create_domino, tiny_mod_args),
        "Excluding 3 interactions that do not pair one receptor \\(R\\) with one ligand \\(L\\): L L, R R, NA L"
    )
    # the excluded rows do not add receptors or linkages; the result matches the valid rows alone
    expect_equal(dom@linkages, tiny_created_dom1@linkages)
    expect_equal(dom@cor, tiny_created_dom1@cor)

    no_rl <- rl_map_tiny
    no_rl$type_A <- "L"
    no_rl$type_B <- "L"
    tiny_mod_args$rl_map <- no_rl
    expect_error(
        suppressWarnings(do.call(create_domino, tiny_mod_args)),
        "No rl_map rows pair a receptor \\(R\\) with a ligand \\(L\\)"
    )
})
