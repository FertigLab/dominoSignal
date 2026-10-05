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
    # stored fixture predates version recording, so compare without it
    expect_equal(dom@misc$create_version, as.character(packageVersion("dominoSignal")))
    dom@misc$create_version <- NULL
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
            min_tf_pval = 0.05,
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
