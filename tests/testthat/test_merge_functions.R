# Write unit tests for each of the merge functions as well as the overall merge function
# Use test objects pbmc_dom_built_tiny and pbmc_dom_built_tiny_2

dom1 <- domino2:::pbmc_dom_built_tiny
dom2 <- domino2:::pbmc_dom_built_tiny_2
doms <- list(dom1, dom2)
misc_list <- lapply(doms, slot, "misc")
link_list <- lapply(doms, slot, "linkages")

test_that("merge db info works", {
    expect_equal(merge_db_info(doms), dom1@db_info)
})

test_that("merge_misc_build_create works", {
    expect_equal(merge_misc_build_create(misc_list), list("build" = TRUE, "create" = TRUE))
})

test_that("merge_misc_rl_map works", {
    expect_equal(merge_misc_rl_map(misc_list), dom1@misc$rl_map)
})


test_that("merge_misc_cor works", {
    expect_equal(
        merge_misc_cor(misc_list, "mean"),
        structure(c(
            -0.0256397186375477, 0, 0, 0, 0.249066622024575,
            0, 0, 0, -0.0112478347367109, 0, 0, 0
        ), dim = 4:3, dimnames = list(
            c("IL7R", "IL2RG", "FAS", "FCER2"), c("ZNF257", "ATF4", "RUNX1")
        ))
    )
    expect_equal(
        merge_misc_cor(misc_list, "range"),
        structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), dim = 4:3, dimnames = list(
            c("IL7R", "IL2RG", "FAS", "FCER2"), c("ZNF257", "ATF4", "RUNX1")
        ))
    )
})

test_that("merge_misc_percent works", {
    expect_equal(merge_misc_percent(misc_list, "mean"), dom1@misc$cl_rec_percent)
    expect_equal(merge_misc_percent(misc_list, "range"), structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        dim = 4:3, dimnames = list(
            c("IL7R", "IL2RG", "FAS", "FCER2"), c(
                "CD8_T_cell", "CD14_monocyte",
                "B_cell"
            )
        )
    ))
})

test_that("merge_misc_build works", {
    expect_equal(merge_misc_build(misc_list), list(
        max_tf_per_clust = Inf, min_tf_pval = c(0.05, 0.15), max_rec_per_tf = Inf,
        rec_tf_cor_threshold = c(0.1, 0.075), min_rec_percentage = c(
            0.01,
            0.02
        )
    ))
})

test_that("merge_misc works", {
    expect_equal(merge_misc(doms, "range"), list(
        build = TRUE, create = TRUE, rl_map = structure(list(
            R.gene = c(
                "FCER2",
                "IGF1R", "FAS", "IL7R,IL2RG"
            ), L.gene = c(
                "ITGB2,ITGAM", "IGF1",
                "TNF", "IL7"
            ), R.uniprot = c("P06734", "P08069", "P25445", "P16871,P31785"),
            L.uniprot = c("P05107,P11215", "P05019", "P01375", "P13232"),
            R.name = c("FCER2", "IGF1R", "FAS", "IL7_receptor"), L.name = c(
                "integrin_aMb2_complex",
                "IGF1", "TNF", "IL7"
            )
        ), row.names = c(NA, -4L), class = "data.frame"),
        rec_cor = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), dim = 4:3, dimnames = list(
            c("IL7R", "IL2RG", "FAS", "FCER2"), c(
                "ZNF257", "ATF4",
                "RUNX1"
            )
        )), cl_rec_percent = structure(c(
            0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ), dim = 4:3, dimnames = list(c(
            "IL7R",
            "IL2RG", "FAS", "FCER2"
        ), c(
            "CD8_T_cell", "CD14_monocyte",
            "B_cell"
        ))), build_vars = list(max_tf_per_clust = Inf, min_tf_pval = c(
            0.05,
            0.15
        ), max_rec_per_tf = Inf, rec_tf_cor_threshold = c(
            0.1,
            0.075
        ), min_rec_percentage = c(0.01, 0.02))
    ))
    expect_equal(merge_misc(doms, "mean"), list(
        build = TRUE, create = TRUE, rl_map = structure(list(
            R.gene = c(
                "FCER2",
                "IGF1R", "FAS", "IL7R,IL2RG"
            ), L.gene = c(
                "ITGB2,ITGAM", "IGF1",
                "TNF", "IL7"
            ), R.uniprot = c("P06734", "P08069", "P25445", "P16871,P31785"),
            L.uniprot = c("P05107,P11215", "P05019", "P01375", "P13232"),
            R.name = c("FCER2", "IGF1R", "FAS", "IL7_receptor"), L.name = c(
                "integrin_aMb2_complex",
                "IGF1", "TNF", "IL7"
            )
        ), row.names = c(NA, -4L), class = "data.frame"),
        rec_cor = structure(c(
            -0.0256397186375477, 0, 0, 0, 0.249066622024575,
            0, 0, 0, -0.0112478347367109, 0, 0, 0
        ), dim = 4:3, dimnames = list(
            c("IL7R", "IL2RG", "FAS", "FCER2"), c(
                "ZNF257", "ATF4",
                "RUNX1"
            )
        )), cl_rec_percent = structure(c(
            0.216666666666667,
            0.416666666666667, 0.0333333333333333, 0.116666666666667,
            0.291666666666667, 0.358333333333333, 0.0166666666666667,
            0.125, 0.241666666666667, 0.308333333333333, 0.0416666666666667,
            0.158333333333333
        ), dim = 4:3, dimnames = list(c(
            "IL7R",
            "IL2RG", "FAS", "FCER2"
        ), c(
            "CD8_T_cell", "CD14_monocyte",
            "B_cell"
        ))), build_vars = list(max_tf_per_clust = Inf, min_tf_pval = c(
            0.05,
            0.15
        ), max_rec_per_tf = Inf, rec_tf_cor_threshold = c(
            0.1,
            0.075
        ), min_rec_percentage = c(0.01, 0.02))
    ))
})

test_that("merge_counts works", {
    expect_equal(
        rownames(merge_counts(doms, "intersect")),
        c(
            "IL7R", "TNF", "IL2RG", "IL7", "FAS", "IGF1R", "ITGAM",
            "FCER2", "ITGB2"
        )
    )
    expect_null(rownames(merge_counts(doms, "outersect")))
})

test_that("merge_de works", {
    expect_equal(merge_de(doms, "mean"), dom1@clust_de)
    expect_equal(merge_de(doms, "range"), structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0), dim = c(3L, 3L), dimnames = list(
        c("ZNF257", "ATF4", "RUNX1"), c(
            "CD8_T_cell", "CD14_monocyte",
            "B_cell"
        )
    )))
})

test_that("merge_zs works", {
    expect_null(rownames(merge_zs(doms, "outersect")))
    expect_equal(rownames(merge_zs(doms, "intersect")), c(
        "IL7R", "TNF", "IL2RG", "IL7", "FAS", "IGF1R", "ITGAM",
        "FCER2", "ITGB2"
    ))
})

test_that("merge_cor works", {
    expect_equal(merge_cor(doms, "mean"), structure(c(
        0.120194321460811, 0, 0.0622076176469977, -0.00489304790873231,
        0.080524889725313, 0, 0.0957948833946582, 0.217338454759885,
        0.0518524925657955, 0, 0.108889086402934, 0.00664182136694283
    ), dim = 4:3, dimnames = list(c("FCER2", "IGF1R", "FAS", "IL7_receptor"), c("ZNF257", "ATF4", "RUNX1"))))
    expect_equal(merge_cor(doms, "range"), structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), dim = 4:3, dimnames = list(
        c("FCER2", "IGF1R", "FAS", "IL7_receptor"), c(
            "ZNF257", "ATF4",
            "RUNX1"
        )
    )))
})

test_that("merge_clusters works", {
    expect_equal(class(merge_clusters(doms)), "factor")
    expect_equal(length(merge_clusters(doms)), sum(sapply(doms, function(x) length(x@clusters))))
})

test_that("merge_features works", {
    expect_equal(dim(merge_features(doms, "outersect")), c(0, sum(sapply(doms, function(x) ncol(x@features)))))
    expect_equal(dim(merge_features(doms, "intersect")), c(3, sum(sapply(doms, function(x) ncol(x@features)))))
})

test_that("merge_link_complex works", {
    expect_equal(merge_link_complex(link_list), dom1@linkages$complex)
})

test_that("merge_link_rl works", {
    expect_equal(merge_link_rl(link_list), dom1@linkages$rec_lig)
})

test_that("merge_link_regulon works", {
    expect_equal(names(merge_link_regulon(link_list)), c("ATF4", "ZNF257", "RUNX1"))
    expect_equal(class(merge_link_regulon(link_list)), "list")
})

test_that("merge_link_tfr works", {
    expect_equal(merge_link_tfr(link_list), list(ZNF257 = "FCER2", ATF4 = c("IL7_receptor", "FAS", "FCER2"), RUNX1 = "FAS"))
})

test_that("merge_link_in_lig works", {
    expect_equal(merge_link_in_lig(link_list), list(CD8_T_cell = "integrin_aMb2_complex", CD14_monocyte = c(
        "IL7",
        "integrin_aMb2_complex"
    ), B_cell = "TNF"))
})

test_that("merge_link_clust_tf works", {
    expect_equal(merge_link_clust_tf(link_list, "union"), list(CD8_T_cell = "ZNF257", CD14_monocyte = "ATF4", B_cell = "RUNX1"))
    expect_equal(merge_link_clust_tf(link_list, "outersect"), list(
        CD8_T_cell = character(0), CD14_monocyte = character(0),
        B_cell = character(0)
    ))
})

test_that("merge_link_clust_tf_rec works", {
    expect_equal(merge_link_clust_tf_rec(link_list, "union"), list(CD8_T_cell = list(ZNF257 = "FCER2"), CD14_monocyte = list(
        ATF4 = "IL7_receptor"
    ), B_cell = list(RUNX1 = "FAS")))
    expect_equal(merge_link_clust_tf_rec(link_list, "outersect"), list(CD8_T_cell = list(), CD14_monocyte = list(), B_cell = list()))
})

test_that("merge_link_clust_rec works", {
    expect_equal(merge_link_clust_rec(link_list, "union"), list(CD8_T_cell = "FCER2", CD14_monocyte = c(
        "IL7_receptor",
        "FCER2"
    ), B_cell = "FAS"))
    expect_equal(merge_link_clust_rec(link_list, "outersect"), list(
        CD8_T_cell = character(0), CD14_monocyte = "FCER2",
        B_cell = character(0)
    ))
})

test_that("merge_linkages works", {
    expect_equal(names(merge_linkages(doms, "union")), c(
        "complexes", "rec_lig", "tf_targets", "clust_tf", "tf_rec",
        "clust_tf_rec", "clust_rec", "clust_incoming_lig"
    ))
    expect_equal(class(merge_linkages(doms, "union")), "list")
})

test_that("merge_cl_signaling works", {
    expect_equal(merge_cl_signaling(doms, "mean"), list(
        CD8_T_cell = structure(c(
            0, 0.282751995942753, 0, 0.1520266905147,
            0, 0.0125341723469197
        ), dim = 2:3, dimnames = list(c(
            "ITGB2",
            "ITGAM"
        ), c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell"))),
        CD14_monocyte = structure(c(
            0.235814498860943, 0.0706879989856883,
            0, 0.038006672628675, 0.0854315307031852, 0.00313354308672992
        ), dim = 2:3, dimnames = list(
            c("IL7", "integrin_aMb2_complex"),
            c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell")
        )), B_cell = structure(c(
            0.083994790240008,
            0, 0
        ), dim = c(1L, 3L), dimnames = list("TNF", c(
            "L_CD8_T_cell",
            "L_CD14_monocyte", "L_B_cell"
        )))
    ))
    expect_equal(merge_cl_signaling(doms, "range"), list(
        CD8_T_cell = structure(c(
            0, 0.282751995942753, 0, 0.1520266905147,
            0, 0.0125341723469197
        ), dim = 2:3, dimnames = list(c(
            "ITGB2",
            "ITGAM"
        ), c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell"))),
        CD14_monocyte = structure(c(
            0.235814498860943, 0.141375997971377,
            0, 0.07601334525735, 0.0854315307031852, 0.00626708617345984
        ), dim = 2:3, dimnames = list(
            c("IL7", "integrin_aMb2_complex"),
            c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell")
        )), B_cell = structure(c(
            0.083994790240008,
            0, 0
        ), dim = c(1L, 3L), dimnames = list("TNF", c(
            "L_CD8_T_cell",
            "L_CD14_monocyte", "L_B_cell"
        )))
    ))
})

test_that("merge_signaling works", {
    expect_equal(merge_signaling(doms, "mean"), structure(c(
        0.282751995942753, 0.18859524841616, 0, 0.1520266905147,
        0.038006672628675, 0, 0.0125341723469197, 0.0458493084383225,
        0
    ), dim = c(3L, 3L), dimnames = list(c(
        "R_CD8_T_cell", "R_CD14_monocyte",
        "R_B_cell"
    ), c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell"))))
    expect_equal(merge_signaling(doms, "range"), structure(c(
        0.282751995942753, 0.377190496832319, 0, 0.1520266905147,
        0.07601334525735, 0, 0.0125341723469197, 0.091698616876645, 0
    ), dim = c(3L, 3L), dimnames = list(c(
        "R_CD8_T_cell", "R_CD14_monocyte",
        "R_B_cell"
    ), c("L_CD8_T_cell", "L_CD14_monocyte", "L_B_cell"))))
})

test_that("merge_dom works", {
    expect_equal(class(merge_dom(doms, "intersect", "mean")), structure("domino", package = "domino2"))
    expect_equal(
        slotNames(merge_dom(doms, "intersect", "mean")),
        c(
            "db_info", "z_scores", "counts", "clusters", "features", "cor",
            "linkages", "clust_de", "misc", "cl_signaling_matrices", "signaling"
        )
    )
})