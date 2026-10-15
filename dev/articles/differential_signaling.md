# Differential Signaling Workflow

In this tutorial, we will walk through the workflow for performing
differential signaling analysis using the dominoSignal package. While
the other tutorials focus on signaling networks inferred from a single
data set, many studies involve multiple subjects or conditions,
requiring the comparison of signaling across these groups. The
dominoSignal package implements the Differential Cell Signaling Test
(DCST), a framework for identifying receptor, ligand, and transcription
factor (TF) linkages that are significantly differentially present
across groups of subjects. This method is described in detail in
[Differential cell signaling testing for cell-cell communication
inference from single-cell data by
dominoSignal](https://doi.org/10.1093/bioinformatics/btag089) (2026)
published in Bioinformatics. The publication also includes simulation
benchmarking and application of the method to pancreatic cancer and
immunotherapy data sets. Refer to the publication for full details on
the method and biological validation of the workflow.

## Workflow overview

The differential signaling workflow consists of the following steps:

1.  Build a
    [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
    object per subject to be compared (see the [Getting
    Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html)
    tutorial for details on building a
    [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
    object).
2.  Compile linkages across the subjects into a [linkage_summary
    object](https://fertiglab.github.io/dominoSignal/reference/linkage_summary-class.html)
    using the
    [`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)
    function.
3.  \[Optional\] Subset the `linkage_summary` object to include only the
    subjects of interest.
4.  Create a table of linkage observations within groups of subjects
    using the
    [`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
    function.
5.  Test differences in linkage frequency across groups using the
    [`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
    function, which implements a Fisher’s exact test.
6.  Visualize the results using the
    [`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md)
    function.

The DCST workflow can be performed with 2 or more groups of subjects. In
this tutorial, we will demonstrate the workflow for both a two-group
comparison (condition A vs condition B) and a comparison across more
than two groups.

## Demonstration data

In this tutorial, we use the
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object built in the [Getting
Started](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html)
tutorial. If you have not yet built a
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object, you can do so by following the instructions on that page or use
the code below to load the pre-built
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object.

Instructions to load data

``` r

set.seed(42)
library(dominoSignal)
library(dplyr)
library(kableExtra)

# BiocFileCache helps with file management across sessions
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
data_url <- "https://zenodo.org/records/10951634/files/pbmc_domino_built.rds"
tmp_path <- BiocFileCache::bfcrpath(bfc, data_url)
dom <- readRDS(tmp_path)
```

Because the differential signaling workflow compares linkages across
multiple
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects (one per subject), we will create six synthetic
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects from the same PBMC dataset. For each synthetic subject, we
randomly retain a subset of the cells without replacement while keeping
their true cluster labels and then perturb the cluster labels of a
subset of those cells by reshuffling within that subset. This approach
preserves some amount of per-cluster signaling while making sure each
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object is distinct from the others, giving the example workflow
something to compare.[^1]

Note that reshuffling is only for demonstration purposes and does not
represent any biological differences. In practice, you would build one
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object from each biological subject or sample to be used in the
comparison. To ensure results are comparable across subjects, use the
same parameters in
[`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
and
[`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
for all subjects. Refer to the [DCST
publication](https://doi.org/10.1093/bioinformatics/btag089) for an
example of this workflow applied to real data.

Creating synthetic replicates

To start, we extract the inputs from the original
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object so we can supply them to the
[`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
function for the synthetic subjects.

``` r

rl_map <- dom_database(dom, name_only = FALSE)
features <- dom_tf_activation(dom)
# Access counts slot directly to maintain sparse matrix format
counts <- slot(dom, "counts")
zscores <- dom_zscores(dom)
tf_targets <- dom_linkages(dom, link_type = "tf-target")
clusters <- dom_clusters(dom, labels = TRUE)
```

We can use a helper function to resample the cells and shuffle the
cluster labels for a subset of them, and then rebuild the
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
object with
[`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
and
[`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md)
using the same parameters across all subjects so that they can be
compared.[^2]

``` r

build_reps <- function(clusters, seed, keep_frac, shuffle_frac) {
    # Local seed for reproducibility without affecting global seed
    withr::with_seed(seed, {
        n_cells <- length(clusters)
        keep_cells <- sort(sample.int(n_cells, size = round(n_cells * keep_frac),
            replace = FALSE))

        resampled_clusters <- clusters[keep_cells]
        features <- features[, keep_cells, drop = FALSE]
        counts <- counts[, keep_cells, drop = FALSE]
        zscores <- zscores[, keep_cells, drop = FALSE]
        shuffled_idx <- sample(seq_along(resampled_clusters), size = round(shuffle_frac *
            length(resampled_clusters)), replace = FALSE)

        resampled_clusters[shuffled_idx] <- sample(resampled_clusters[shuffled_idx])

        dom_shuffle <- create_domino(rl_map = rl_map, features = features, counts = counts,
            z_scores = zscores, clusters = resampled_clusters, tf_targets = tf_targets,
            use_complexes = TRUE, remove_rec_dropout = FALSE, verbose = FALSE)

        dom_shuffle <- build_domino(dom = dom_shuffle, max_tf_pval = 0.05, max_tf_per_clust = 25,
            max_rec_per_tf = 25, rec_tf_cor_threshold = 0.1, min_rec_percentage = 0.05)
    })

    return(dom_shuffle)
}

params <- data.frame(set = c(161803, 314159, 6631034, 6021023, 831446, 98067), keep = c(0.95,
    0.95, 0.95, 0.95, 0.9, 0.9), shuffle = c(0.15, 0.1, 0.15, 0.1, 0.15, 0.1))

# This does take a few minutes to run
dom_list <- purrr::pmap(params, function(set, keep, shuffle) {
    build_reps(clusters, seed = set, keep_frac = keep, shuffle_frac = shuffle)
})

names(dom_list) <- paste0("P", seq_along(dom_list))
```

We now have six
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects that we can treat as different subjects (P1-6).

## The `linkage_summary` class

The
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
class is the structure used for the differential signaling workflow. It
contains linkages present in each cluster across the
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects it contains, along with subject metadata. The object contains
three slots:

- `subject_names`: a factor with a unique name for each subject
  ([`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object) represented in the summary.
- `subject_meta`: a data frame with one row per subject including
  provided metadata (such as condition, treatment, or other covariates)
  that can be used to group or subset subjects for comparison.
- `subject_linkages`: a nested list organized as subject to cluster to
  linkage type (`tfs` for active transcription factors, `rec` for active
  receptors, `incoming_lig` for ligands that can activate the cluster,
  `tfs_rec` for transcription factor to receptor pairs, and `rec_lig`
  for receptor to ligand pairs).

To create a
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
object, we use the
[`summarize_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/summarize_linkages.md)
function, which takes a named list of
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects and a data frame of subject metadata. The first column of the
data frame must contain subject names that match the names of the
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects in the list.

``` r

subject_metas <- data.frame(patient = c("P1", "P2", "P3", "P4", "P5", "P6"), time = c("pre",
    "post", "pre", "post", "pre", "post"), condition = c("control", "treatment",
    "treatment", "control", "control", "control"), batch = c("A", "A", "B", "B",
    "C", "C"), stringsAsFactors = FALSE)

link_summary <- summarize_linkages(dom_list, subject_metas, subject_names = subject_metas$patient)

print(link_summary)
#> A linkage summary object of 6 subjects with 4 metadata annotations and linkages between 9 clusters.
```

The `print` and `show` methods for the
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
class report the number of subjects, the number of metadata annotations,
and the number of clusters present in the summary.

Since we built this linkage summary from six
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects with six distinct groups and two distinct conditions, we can
subset the summary to include only the subjects of interest.

To create a comparison for only the treatment condition, we can subset
the
[`linkage_summary()`](https://FertigLab.github.io/dominoSignal/dev/reference/linkage_summary-class.md)
object by the relevant metadata column.

``` r

link_treat <- subset(link_summary, subset = condition == "treatment")
print(link_treat)
#> A linkage summary object of 2 subjects with 4 metadata annotations and linkages between 9 clusters.
```

Subjects can also be selected by name or with compound expressions.

``` r

subset(link_summary, subset = subject_names %in% c("P1", "P3"))
#> A linkage summary object of 2 subjects with 4 metadata annotations and linkages between 9 clusters.
subset(link_summary, subset = condition == "treatment" | time == "pre")
#> A linkage summary object of 4 subjects with 4 metadata annotations and linkages between 9 clusters.
```

As an alternative to creating separate subset objects, the other
functions in the differential signaling workflow accept a
`subject_names` argument that can be used to restrict the subjects
considered in the analysis without needing to create a separate
`linkage_summary` object.

## Counting linkages across groups

The
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
function tabulates the number of subjects with a given linkage in a
specified cluster, which can also be broken down by a metadata column.
The `linkage` argument selects which of the linkage types to count.

### Total counts

If no grouping variable is specified and no subject names are provided,
the function will return a data frame with the total number of subjects
with each linkage in the selected cluster.

``` r

count_all <- count_linkage(link_summary, cluster = "CD8_T_cell", linkage = "rec_lig") %>%
    slice_sample(n = 5)
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(count_all, caption = "Total counts of receptor-ligand linkages in CD8 T Cells"),
    bootstrap_options = c("striped", "hover"))
```

| feature          | total_count |
|:-----------------|------------:|
| TNFRSF1B \<- LTA |           1 |
| CD2 \<- CD58     |           6 |
| LAMP1 \<- FAM3C  |           6 |
| CXCR3 \<- CXCL11 |           6 |
| TNFRSF1B \<- GRN |           1 |

Total counts of receptor-ligand linkages in CD8 T Cells {.table .table
.table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

If specific subjects are of interest, the counts can be restricted to
those subjects by providing a vector of the subject names.

``` r

count_subjects <- count_linkage(link_summary, cluster = "CD8_T_cell", linkage = "tfs",
    subject_names = c("P1", "P2"))
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(count_subjects[1:5, ], caption = "Subset of counts of TF linkages in patients P1 & P2 in CD8 T Cells"),
    bootstrap_options = c("striped", "hover"))
```

| feature | total_count |
|:--------|------------:|
| ATF6    |           2 |
| CHURC1  |           1 |
| CREM    |           2 |
| DBP     |           2 |
| ELF3    |           2 |

Subset of counts of TF linkages in patients P1 & P2 in CD8 T Cells
{.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

### Grouped counts

The
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
function can also return counts that are grouped by a metadata column.
For example, we can count the linkages for subjects in the treatment
condition grouped by the “time” column to see how many subjects in each
time point have each linkage.

``` r

counts_time_trt <- count_linkage(link_treat, cluster = "CD14_monocyte", linkage = "tfs_rec",
    group.by = "time")
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(slice_sample(counts_time_trt, n = 5), caption = "Counts of TF-Receptor linkages in CD14 Monocytes grouped by time for treatment subjects"),
    bootstrap_options = c("striped", "hover"))
```

| feature            | total_count | post | pre |
|:-------------------|------------:|-----:|----:|
| CEBPE \<- TNFRSF1A |           2 |    1 |   1 |
| CEBPA \<- CSF1R    |           2 |    1 |   1 |
| CEBPE \<- CSF1R    |           2 |    1 |   1 |
| IRF8 \<- TNFRSF1B  |           2 |    1 |   1 |
| ATF3 \<- TNFRSF1A  |           2 |    1 |   1 |

Counts of TF-Receptor linkages in CD14 Monocytes grouped by time for
treatment subjects {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

This is equivalent to grouping by the “time” column and providing the
subject names of only the subjects in the treatment condition.

``` r

counts_time_subj <- count_linkage(link_summary, cluster = "CD14_monocyte", linkage = "tfs_rec",
    group.by = "time", subject_names = c("P2", "P3"))
identical(counts_time_trt, counts_time_subj)
#> [1] TRUE
```

Of course, more than one subject can be in each group, and the counts
will reflect the number of subjects in each group with each linkage.

``` r

counts_condition <- count_linkage(link_summary, cluster = "naive_CD4_T_cell", linkage = "incoming_lig",
    group.by = "condition") %>%
    slice_sample(n = 5)
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(counts_condition, caption = "Counts of incoming ligand linkages in naive CD4 T Cells grouped by condition"),
    bootstrap_options = c("striped", "hover"))
```

| feature | total_count | control | treatment |
|:--------|------------:|--------:|----------:|
| DHCR24  |           6 |       4 |         2 |
| CD70    |           6 |       4 |         2 |
| CDH1    |           6 |       4 |         2 |
| IL7     |           6 |       4 |         2 |
| PTGES2  |           1 |       1 |         0 |

Counts of incoming ligand linkages in naive CD4 T Cells grouped by
condition {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

And more than two groups can be compared at once.

``` r

counts_batch <- count_linkage(link_summary, cluster = "CD14_monocyte", linkage = "rec",
    group.by = "batch") %>%
    slice_sample(n = 5)
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(counts_batch, caption = "Counts of receptor linkages in CD14 Monocytes grouped by batch"),
    bootstrap_options = c("striped", "hover"))
```

| feature | total_count |   A |   B |   C |
|:--------|------------:|----:|----:|----:|
| MARCO   |           6 |   2 |   2 |   2 |
| C5AR1   |           6 |   2 |   2 |   2 |
| LTBR    |           4 |   2 |   0 |   2 |
| NOTCH2  |           6 |   2 |   2 |   2 |
| CD36    |           6 |   2 |   2 |   2 |

Counts of receptor linkages in CD14 Monocytes grouped by batch {.table
.table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

## Testing differential linkages

The
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
function tabulates the number of subjects with a given linkage in a
specified cluster, but it does not test whether the frequency of that
linkage is significantly different across groups. The
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
function uses a contingency table of presence/absence of each linkage
across groups to perform a Fisher’s exact test. The p-values across all
tested linkages are adjusted for multiple comparisons using the false
discovery rate method.

The
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
function returns a data frame that includes the cluster, linkage,
grouping variable, test name (always fisher’s exact test), the linkage
being compared, the test statistics (odds ratio, if applicable, and
p-value, both raw and adjusted), total number of subjects where linkage
is active, number of subjects in each group where linkage is active,
number of total subjects compared, and the total number of subjects in
each group.

### Two-group comparisons

When exactly two groups are compared, the
[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) function
returns an odds ratio in addition to the p-value, reflecting the
relative odds of the linkage being present in one group compared to the
other.

``` r

dcst_time_trt <- test_differential_linkages(link_treat, cluster = "CD14_monocyte",
    linkage = "tfs_rec", group.by = "time")
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(slice_sample(dcst_time_trt, n = 5), caption = "DCST of TF-Receptor linkages in CD14 Monocytes grouped by time for treatment subjects"),
    bootstrap_options = c("striped", "hover"))
```

|  | cluster | linkage | group.by | test_name | feature | odds.ratio | p.value | p.adj | total_count | post_count | pre_count | total_n | post_n | pre_n |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ZNF77 \<- CD36 | CD14_monocyte | tfs_rec | time | fishers.exact | ZNF77 \<- CD36 | 0 | 1 | 1 | 2 | 1 | 1 | 2 | 1 | 1 |
| FOSB \<- CD86 | CD14_monocyte | tfs_rec | time | fishers.exact | FOSB \<- CD86 | 0 | 1 | 1 | 2 | 1 | 1 | 2 | 1 | 1 |
| FOSL1 \<- TNFRSF1B | CD14_monocyte | tfs_rec | time | fishers.exact | FOSL1 \<- TNFRSF1B | 0 | 1 | 1 | 2 | 1 | 1 | 2 | 1 | 1 |
| KLF4 \<- FcRn_complex | CD14_monocyte | tfs_rec | time | fishers.exact | KLF4 \<- FcRn_complex | Inf | 1 | 1 | 1 | 1 | 0 | 2 | 1 | 1 |
| ZNF516 \<- FCGR2A | CD14_monocyte | tfs_rec | time | fishers.exact | ZNF516 \<- FCGR2A | 0 | 1 | 1 | 2 | 1 | 1 | 2 | 1 | 1 |

DCST of TF-Receptor linkages in CD14 Monocytes grouped by time for
treatment subjects {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

As with the
[`count_linkage()`](https://FertigLab.github.io/dominoSignal/dev/reference/count_linkage.md)
function, the same results can be obtained by providing the subject
names of only the subjects in the treatment condition.

``` r

dcst_time_subj <- test_differential_linkages(link_summary, cluster = "CD14_monocyte",
    linkage = "tfs_rec", group.by = "time", subject_names = c("P2", "P3"))
identical(dcst_time_trt, dcst_time_subj)
#> [1] TRUE
```

The
[`test_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/test_differential_linkages.md)
function does require that any grouping variable has at least two groups
with at least one subject in each group.

``` r

test_differential_linkages(link_treat, cluster = "CD14_monocyte", linkage = "tfs_rec",
    group.by = "condition")
#> Error in `test_differential_linkages()`:
#> ! At least 2 groups of condition must be represented among the included subjects; found: treatment
```

It is obviously preferable to have more than one subject in each group
(though the function does return results for groups with only one
subject). Note that there are no significant differences between the
groups in this example, which is expected since the subjects were
generated from the same data set.

``` r

dcst_time <- test_differential_linkages(link_summary, cluster = "CD8_T_cell", linkage = "rec",
    group.by = "time") %>%
    dplyr::arrange(desc(odds.ratio))
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(slice_sample(dcst_time, n = 5), caption = "DCST of receptor-ligand linkages in CD8 T Cells grouped by time"),
    bootstrap_options = c("striped", "hover"))
```

|  | cluster | linkage | group.by | test_name | feature | odds.ratio | p.value | p.adj | total_count | post_count | pre_count | total_n | post_n | pre_n |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CD27 | CD8_T_cell | rec | time | fishers.exact | CD27 | 0 | 1 | 1 | 6 | 3 | 3 | 6 | 3 | 3 |
| SORL1 | CD8_T_cell | rec | time | fishers.exact | SORL1 | 0 | 1 | 1 | 6 | 3 | 3 | 6 | 3 | 3 |
| CD8_receptor | CD8_T_cell | rec | time | fishers.exact | CD8_receptor | 0 | 1 | 1 | 6 | 3 | 3 | 6 | 3 | 3 |
| CXCR4 | CD8_T_cell | rec | time | fishers.exact | CXCR4 | 0 | 1 | 1 | 6 | 3 | 3 | 6 | 3 | 3 |
| PTGDR | CD8_T_cell | rec | time | fishers.exact | PTGDR | 0 | 1 | 1 | 6 | 3 | 3 | 6 | 3 | 3 |

DCST of receptor-ligand linkages in CD8 T Cells grouped by time {.table
.table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

### Multi-group comparisons

When more than two groups are compared, the
[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) function
does not return an odds ratio (so it is set to `NA` in the DCST output),
but the p-value still reflects whether the frequency of the linkage is
significantly different across groups.

``` r

dcst_batch <- test_differential_linkages(link_summary, cluster = "naive_CD4_T_cell",
    linkage = "rec", group.by = "batch") %>%
    arrange(p.value)
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(slice_sample(dcst_batch, n = 5), caption = "DCST of receptors in naive CD4 T Cells grouped by batch"),
    bootstrap_options = c("striped", "hover"))
```

|  | cluster | linkage | group.by | test_name | feature | odds.ratio | p.value | p.adj | total_count | A_count | B_count | C_count | total_n | A_n | B_n | C_n |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CD44 | naive_CD4_T_cell | rec | batch | fishers.exact | CD44 | NA | 1 | 1 | 6 | 2 | 2 | 2 | 6 | 2 | 2 | 2 |
| PTGER2 | naive_CD4_T_cell | rec | batch | fishers.exact | PTGER2 | NA | 1 | 1 | 1 | 0 | 0 | 1 | 6 | 2 | 2 | 2 |
| TNFRSF25 | naive_CD4_T_cell | rec | batch | fishers.exact | TNFRSF25 | NA | 1 | 1 | 6 | 2 | 2 | 2 | 6 | 2 | 2 | 2 |
| RORA | naive_CD4_T_cell | rec | batch | fishers.exact | RORA | NA | 1 | 1 | 6 | 2 | 2 | 2 | 6 | 2 | 2 | 2 |
| HLA-C | naive_CD4_T_cell | rec | batch | fishers.exact | HLA-C | NA | 1 | 1 | 6 | 2 | 2 | 2 | 6 | 2 | 2 | 2 |

DCST of receptors in naive CD4 T Cells grouped by batch {.table .table
.table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

There is no limit to the number of groups that can be compared at once
(as long as there is at least one subject in each group).

``` r

dcst_patient <- test_differential_linkages(link_summary, cluster = "naive_CD4_T_cell",
    linkage = "rec_lig", group.by = "patient") %>%
    arrange(p.value)
```

Show table code

``` r

kableExtra::kable_styling(knitr::kable(slice_sample(dcst_patient, n = 5), caption = "DCST of receptor-ligand linkages in naive CD4 T Cells grouped by patient"),
    bootstrap_options = c("striped", "hover"))
```

|  | cluster | linkage | group.by | test_name | feature | odds.ratio | p.value | p.adj | total_count | P1_count | P2_count | P3_count | P4_count | P5_count | P6_count | total_n | P1_n | P2_n | P3_n | P4_n | P5_n | P6_n |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IL27_receptor \<- IL27 | naive_CD4_T_cell | rec_lig | patient | fishers.exact | IL27_receptor \<- IL27 | NA | 1 | 1 | 4 | 1 | 1 | 1 | 1 | 0 | 0 | 6 | 1 | 1 | 1 | 1 | 1 | 1 |
| NCR3 \<- NCR3LG1 | naive_CD4_T_cell | rec_lig | patient | fishers.exact | NCR3 \<- NCR3LG1 | NA | 1 | 1 | 5 | 1 | 1 | 1 | 0 | 1 | 1 | 6 | 1 | 1 | 1 | 1 | 1 | 1 |
| ICOS \<- TNF | naive_CD4_T_cell | rec_lig | patient | fishers.exact | ICOS \<- TNF | NA | 1 | 1 | 5 | 0 | 1 | 1 | 1 | 1 | 1 | 6 | 1 | 1 | 1 | 1 | 1 | 1 |
| LAMP1 \<- FAM3C | naive_CD4_T_cell | rec_lig | patient | fishers.exact | LAMP1 \<- FAM3C | NA | 1 | 1 | 6 | 1 | 1 | 1 | 1 | 1 | 1 | 6 | 1 | 1 | 1 | 1 | 1 | 1 |
| NCR3 \<- BAG6 | naive_CD4_T_cell | rec_lig | patient | fishers.exact | NCR3 \<- BAG6 | NA | 1 | 1 | 5 | 1 | 1 | 1 | 0 | 1 | 1 | 6 | 1 | 1 | 1 | 1 | 1 | 1 |

DCST of receptor-ligand linkages in naive CD4 T Cells grouped by patient
{.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;"}

The [DCST paper](https://doi.org/10.1093/bioinformatics/btag089)
includes additional information on the recommended number of subjects
for robust results (ideally at least 15 per group), as well as an
approach to simulate biological replicates using bootstrapping of cells
for experiments with small numbers of replicates (scripts for the
manuscript can be found at [the corresponding GitHub
repository](https://github.com/FertigLab/Differential_Cell_Signaling_Test)).

## Visualizing differential linkages

In order to visualize the results of the DCST, the
[`plot_differential_linkages()`](https://FertigLab.github.io/dominoSignal/dev/reference/plot_differential_linkages.md)
function can be used to create a heatmap of the test statistic along
with row-level bar plots showing the proportions of subjects in which
linkages are active.

The function has arguments for the test statistic (the column of the
DCST result that contains the statistic to be plotted), the range of
values to be plotted, the heatmap gradient, and the color palette to use
for the groups being compared.

### Two-group comparisons

An example of a two-group comparison is shown below, where the p-value
is used as the test statistic and rows are ordered by ascending p-value.
All p-values are 1 in this example, so every cell shows the `maximum`
gradient color.

``` r

plot_differential_linkages(differential_linkages = slice_sample(dcst_time_trt, n = 15),
    test_statistic = "p.value", stat_range = c(0, 1), stat_ranking = "ascending")
```

![Heatmap of p-values for receptor-ligand linkages in CD8 T cell cluster
differing between pre and post timepoint for treatment condition with
row barplots showing proportion of subjects in each group with active
linkage.](differential_signaling_files/figure-html/plot%20time%20treat-1.png)

``` r

plot_differential_linkages(differential_linkages = slice_sample(dcst_time, n = 15),
    test_statistic = "p.value", gradient = c(minimum = "#c44f4f", maximum = "#ffd771"),
    stat_ranking = "ascending")
```

![Heatmap of p-values for receptors in CD8 T cells differing between pre
and post timepoint with row barplots showing proportion of subjects in
each group with active
receptor.](differential_signaling_files/figure-html/plot%20time-1.png)

If comparing between two groups, it is also possible to use the odds
ratio as the test statistic, which will show the relative odds of a
linkage being present in one group compared to the other.

``` r

plot_differential_linkages(differential_linkages = slice_sample(dcst_time, n = 15),
    test_statistic = "odds.ratio", stat_range = c(0, 4), stat_ranking = "descending",
    gradient = c(minimum = "gray90", maximum = "orange"))
```

![Heatmap of odds ratios for receptors in CD8 T cells differing between
pre and post timepoints with row barplots showing proportion of subjects
in each group with active
receptor.](differential_signaling_files/figure-html/plot%20dcst%20time%20odds-1.png)

### Multi-group comparisons

When comparing between more than two groups, the odds ratio is not
available, so the p-value (or adjusted p-value) is the only option for
the test statistic.

``` r

plot_differential_linkages(differential_linkages = slice_sample(dcst_batch, n = 15),
    test_statistic = "p.value", stat_range = c(0, 1), stat_ranking = "ascending",
    gradient = c(minimum = "#ff4c43", maximum = "#ffe6c0"))
```

![Heatmap of p-values for receptors in naive CD4 T cells differing
between batches with row barplots showing proportion of subjects in each
group with active
receptor.](differential_signaling_files/figure-html/plot%20dcst%20batch-1.png)

A custom color can be assigned to each group’s barplot annotation using
the `group_palette` argument, which takes a named vector of colors where
names correspond to the group names in the DCST result.

``` r

plot_differential_linkages(differential_linkages = slice_sample(dcst_batch, n = 15),
    test_statistic = "p.value", stat_range = c(0, 1), stat_ranking = "ascending",
    gradient = c(minimum = "#7b43ff", maximum = "#c0c5ff"), group_palette = c(A = "#64a5e6",
        B = "#71fcad", C = "#c968c3"))
```

![](differential_signaling_files/figure-html/plot%20dcst%20group%20colors-1.png)

This can be useful for consistency across multiple plots or to match
publication color schemes.

## Example run

Here’s an example of what a workflow might look like from the list of
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects to the final visualization of differential linkages. We begin
with the list of synthetic
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects created above, create the metadata, create a `link_summary()`
object, examine the counts for the linkages,and then run the DCST for a
comparison by batch. In this example, we focus on incoming ligands in
the CD14 monocyte cluster.

``` r

subject_metas <- data.frame(patient = c("P1", "P2", "P3", "P4", "P5", "P6"), time = c("pre",
    "post", "pre", "post", "pre", "post"), condition = c("control", "treatment",
    "treatment", "control", "control", "control"), batch = c("A", "A", "B", "B",
    "C", "C"), stringsAsFactors = FALSE)

link_summary <- summarize_linkages(dom_list, subject_metas, subject_names = subject_metas$patient)

counts_batch <- count_linkage(link_summary, cluster = "CD14_monocyte", linkage = "incoming_lig",
    group.by = "batch")

dcst_batch <- test_differential_linkages(link_summary, cluster = "CD14_monocyte",
    linkage = "incoming_lig", group.by = "batch")

plot_differential_linkages(dcst_batch, test_statistic = "p.value", stat_ranking = "ascending",
    gradient = c(minimum = "#e02c2c", maximum = "#feffcd"), group_palette = c(A = "#2a6aaa",
        B = "#07bd56", C = "#a244a3"))
```

## Related Articles

- [Get started with
  dominoSignal](https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html):
  the basic workflow
- [Interacting with domino
  objects](https://fertiglab.github.io/dominoSignal/articles/domino_object_vignette.html):
  object structure and accessor functions
- [Plotting functions and
  options](https://fertiglab.github.io/dominoSignal/articles/plotting_vignette.html):
  visualization options
- [SCENIC for TF activity
  scoring](https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette.html):
  generating TF activity inputs
- [Using the CellPhoneDB
  database](https://fertiglab.github.io/dominoSignal/articles/cellphonedb_vignette.html):
  preparing a receptor-ligand map

------------------------------------------------------------------------

To report a bug, ask a question, or suggest a feature, please [open an
issue](https://github.com/FertigLab/dominoSignal/issues).

Vignette Build Information

Date last built and session information:

``` r

Sys.Date()
#> [1] "2026-10-15"
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] kableExtra_1.4.1   dplyr_1.2.1        dominoSignal_1.7.1
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.3.0             httr2_1.3.0           formatR_1.14         
#>   [4] biomaRt_2.68.0        rlang_1.3.0           magrittr_2.0.5       
#>   [7] clue_0.3-68           GetoptLong_1.1.1      otel_0.2.0           
#>  [10] matrixStats_1.5.0     compiler_4.6.1        RSQLite_3.53.3       
#>  [13] png_0.1-9             systemfonts_1.3.2     vctrs_0.7.3          
#>  [16] stringr_1.6.0         pkgconfig_2.0.3       shape_1.4.6.1        
#>  [19] crayon_1.5.3          fastmap_1.2.0         backports_1.5.1      
#>  [22] dbplyr_2.6.0          XVector_0.52.0        rmarkdown_2.32       
#>  [25] ragg_1.5.2            purrr_1.2.2           bit_4.6.0            
#>  [28] xfun_0.61             cachem_1.1.0          jsonlite_2.0.0       
#>  [31] progress_1.2.3        blob_1.3.0            broom_1.0.13         
#>  [34] parallel_4.6.1        prettyunits_1.2.0     cluster_2.1.8.2      
#>  [37] R6_2.6.1              bslib_0.12.0          stringi_1.8.9        
#>  [40] RColorBrewer_1.1-3    car_3.1-5             jquerylib_0.1.4      
#>  [43] Rcpp_1.1.2            Seqinfo_1.2.0         iterators_1.0.14     
#>  [46] knitr_1.52            IRanges_2.46.0        Matrix_1.7-5         
#>  [49] igraph_2.3.4          tidyselect_1.2.1      rstudioapi_0.19.0    
#>  [52] abind_1.4-8           yaml_2.3.12           doParallel_1.0.17    
#>  [55] codetools_0.2-20      curl_8.1.0            lattice_0.22-9       
#>  [58] tibble_3.3.1          plyr_1.8.9            withr_3.0.3          
#>  [61] Biobase_2.72.0        KEGGREST_1.52.2       S7_0.2.2             
#>  [64] evaluate_1.0.5        desc_1.4.3            BiocFileCache_3.2.0  
#>  [67] xml2_1.6.0            circlize_0.4.18       Biostrings_2.80.2    
#>  [70] pillar_1.11.1         ggpubr_1.0.0          filelock_1.0.3       
#>  [73] carData_3.0-6         foreach_1.5.2         stats4_4.6.1         
#>  [76] generics_0.1.4        S4Vectors_0.50.3      hms_1.1.4            
#>  [79] ggplot2_4.0.3         scales_1.4.0          glue_1.8.1           
#>  [82] tools_4.6.1           ggsignif_0.6.4        fs_2.1.0             
#>  [85] grid_4.6.1            tidyr_1.3.2           AnnotationDbi_1.74.0 
#>  [88] colorspace_2.1-4      Formula_1.2-6         cli_3.6.6            
#>  [91] textshaping_1.0.5     viridisLite_0.4.3     svglite_2.2.2        
#>  [94] ComplexHeatmap_2.28.0 gtable_0.3.6          rstatix_1.1.0        
#>  [97] sass_0.4.10           digest_0.6.39         BiocGenerics_0.58.1  
#> [100] rjson_0.2.23          htmlwidgets_1.6.4     farver_2.1.2         
#> [103] memoise_2.0.1         htmltools_0.5.9       pkgdown_2.2.1        
#> [106] lifecycle_1.0.5       httr_1.4.9            GlobalOptions_0.1.4  
#> [109] bit64_4.8.6
```

[^1]: Other methods, such as shuffling all cluster labels, results in
    clusters that are all statistically similar which results in few
    linkages reaching significance. While this approach reduces the
    differences between domino objects, it preserves enough signaling to
    infer a reasonable number of (insignificant) linkages.

[^2]: These parameters are fairly permissive to ensure enough linkages
    are inferred for demonstration purposes.
