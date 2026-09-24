test_that("linkage_summary class methods run", {
    expect_s4_class(tiny_linkage_summary, "linkage_summary")
    expect_no_error(print(tiny_linkage_summary))
    expect_no_error(show(tiny_linkage_summary))
    expect_no_error(subset(tiny_linkage_summary, subset = group == "B"))
})

test_that("print reports the object's actual subject, metadata, and cluster counts", {
    expect_output(
        print(tiny_linkage_summary),
        "A linkage summary object of 3 subjects with 2 metadata annotations and linkages between 3 clusters."
    )
})

test_that("show reports the object's actual subject, metadata, and cluster counts", {
    expect_output(
        show(tiny_linkage_summary),
        "A linkage summary object of 3 subjects with 2 metadata annotations and linkages between 3 clusters."
    )
})

test_that("print and show reflect a subsetted object's reduced subject count", {
    sub <- subset(tiny_linkage_summary, subset = group == "A")
    expect_output(
        print(sub),
        "A linkage summary object of 2 subjects with 2 metadata annotations and linkages between 3 clusters."
    )
    expect_output(
        show(sub),
        "A linkage summary object of 2 subjects with 2 metadata annotations and linkages between 3 clusters."
    )
})

test_that("print and show reflect total metadata annotations", {
    big_links <- tiny_linkage_summary
    big_links@subject_meta$extra_col <- c("X", "Y", "Z")
    expect_output(
        print(big_links),
        "A linkage summary object of 3 subjects with 3 metadata annotations and linkages between 3 clusters."
    )
    expect_output(
        show(big_links),
        "A linkage summary object of 3 subjects with 3 metadata annotations and linkages between 3 clusters."
    )
})

test_that("print and show handle empty linkage_summary objects", {
    # prototype object with no slots supplied; S4 skips validity checks in this case
    empty_link <- linkage_summary()
    expect_output(print(empty_link), "An empty linkage summary object \\(0 subjects\\)\\.")
    expect_output(show(empty_link), "An empty linkage summary object \\(0 subjects\\)\\.")
    expect_no_warning(capture.output(show(empty_link)))
})

test_that("show reports 0 clusters rather than -Inf when subjects have no linkages", {
    no_links <- linkage_summary(
        subject_names = factor("dom1"),
        subject_meta = data.frame(ID = "dom1"),
        subject_linkages = list(dom1 = list())
    )
    expect_no_warning(capture.output(show(no_links)))
    expect_output(show(no_links), "linkages between 0 clusters\\.")
})

test_that("linkage_summary validity rejects empty objects", {
    expect_error(
        linkage_summary(
            subject_names = factor(character()),
            subject_meta = data.frame(ID = character()),
            subject_linkages = list()
        ),
        "subject_names must have length > 0"
    )
})

test_that("linkage_summary validity rejects slots of mismatched length", {
    expect_error(
        linkage_summary(
            subject_names = factor(c("dom1", "dom2")),
            subject_meta = data.frame(ID = c("dom1", "dom2")),
            subject_linkages = list(dom1 = list())
        ),
        "subject_linkages must have length 2"
    )
    expect_error(
        linkage_summary(
            subject_names = factor(c("dom1", "dom2")),
            subject_meta = data.frame(ID = "dom1"),
            subject_linkages = list(dom1 = list(), dom2 = list())
        ),
        "subject_meta must have 2 rows"
    )
    expect_error(
        linkage_summary(
            subject_names = factor("dom1"),
            subject_meta = data.frame(row.names = "dom1"),
            subject_linkages = list(dom1 = list())
        ),
        "subject_meta must have at least 1 column"
    )
    expect_true(validObject(tiny_linkage_summary))
})

test_that("subset filters on a subject_meta column", {
    sub <- expect_no_warning(subset(tiny_linkage_summary, subset = group == "A"))
    expect_s4_class(sub, "linkage_summary")
    expect_equal(as.character(sub@subject_names), c("dom1", "dom2"))
    expect_equal(sub@subject_meta$ID, c("dom1", "dom2"))
    expect_named(sub@subject_linkages, c("dom1", "dom2"))
})

test_that("subset filters on subject_names", {
    sub <- subset(tiny_linkage_summary, subset = subject_names %in% c("dom1", "dom3"))
    expect_equal(as.character(sub@subject_names), c("dom1", "dom3"))
    expect_equal(sub@subject_meta$ID, c("dom1", "dom3"))
    expect_named(sub@subject_linkages, c("dom1", "dom3"))
})

test_that("subset filters on a combination of subject_names and subject_meta columns", {
    sub <- subset(tiny_linkage_summary, subset = subject_names %in% c("dom1", "dom2", "dom3") & group == "B")
    expect_equal(as.character(sub@subject_names), "dom3")
    expect_equal(sub@subject_meta$ID, "dom3")
    expect_named(sub@subject_linkages, "dom3")
})

test_that("subset drops unused subject_names factor levels, preserving relative order", {
    sub <- subset(tiny_linkage_summary, subset = group == "A")
    expect_equal(levels(sub@subject_names), c("dom1", "dom2"))
    expect_equal(nlevels(sub@subject_names), length(sub@subject_names))
})

test_that("subset errors when no subjects match", {
    expect_error(
        subset(tiny_linkage_summary, subset = group == "C"),
        "No subjects matched"
    )
    expect_error(
        subset(tiny_linkage_summary, subset = subject_names == "none"),
        "No subjects matched"
    )
})

test_that("subset errors on unknown variables", {
    expect_error(
        subset(tiny_linkage_summary, subset = bogus_var == "A"),
        "unknown variable"
    )
    expect_error(
        subset(tiny_linkage_summary, subset = bogus_var == "A"),
        "bogus_var"
    )
})

test_that("subset does not silently pull unknown variables from the calling environment", {
    # `unrelated` exists in this scope but is not subject_names or a subject_meta column,
    # so it must not be picked up by the subset expression
    unrelated <- c(TRUE, TRUE, FALSE)
    expect_error(
        subset(tiny_linkage_summary, subset = unrelated == TRUE), # nolint
        "unknown variable"
    )
})
