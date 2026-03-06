library(SingleCellExperiment)
library(dominoSignal)

# Note that v0.2.1 is from a previous domino version (v0.2.1); refer to commit d8b0e11 for creation
# Included for regression tests to assess changes to processing functions and outputs
# Assumed to be loaded for saving with other internal data in R/sysdata.rda

# This chunk is identical to the code for exported data
# load data for generation of test results from zenodo repository
# Zenodo host of outputs from SCENIC analysis
data_url <- "https://zenodo.org/records/10951634/files"
temp_dir <- tempdir()

pbmc_dir <- file.path(temp_dir, "/pbmc")
if (!dir.exists(pbmc_dir)) {
    dir.create(pbmc_dir)
}

# SingleCellExperiment object of preprocessed PBMC3K data
download.file(
    url = file.path(data_url, "/pbmc3k_sce.rds"),
    destfile = file.path(pbmc_dir, "/pbmc3k_sce.rds"),
    mode = "wb"
)
pbmc <- readRDS(file.path(pbmc_dir, "/pbmc3k_sce.rds"))

# SCENIC input files
scenic_dir <- file.path(temp_dir, "/scenic")
if (!dir.exists(scenic_dir)) {
    dir.create(scenic_dir)
}
download.file(
    url = file.path(data_url, "/auc_pbmc_3k.csv"),
    destfile = file.path(scenic_dir, "/auc_pbmc_3k.csv"),
    mode = "wb"
)
download.file(
    url = file.path(data_url, "/regulons_pbmc_3k.csv"),
    destfile = file.path(scenic_dir, "/regulons_pbmc_3k.csv"),
    mode = "wb"
)
auc <- read.table(file.path(scenic_dir, "/auc_pbmc_3k.csv"),
    header = TRUE, row.names = 1,
    stringsAsFactors = FALSE, sep = ","
)
regulons <- read.csv(file.path(scenic_dir, "/regulons_pbmc_3k.csv"))

# CellPhoneDB Database
cellphone_url <- "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz"
cellphone_tar <- file.path(temp_dir, "/cellphoneDB_v4.tar.gz")
if (!file.exists(cellphone_tar)) {
    download.file(url = cellphone_url, destfile = cellphone_tar, mode = "wb")
}

cellphone_dir <- file.path(temp_dir, "/cellphoneDB_v4")
untar(tarfile = cellphone_tar, exdir = cellphone_dir)
cellphone_data <- file.path(cellphone_dir, "cellphonedb-data-4.0.0", "data")

interactions <- read.csv(file.path(cellphone_data, "/interaction_input.csv"), stringsAsFactors = FALSE)
complexes <- read.csv(file.path(cellphone_data, "/complex_input.csv"), stringsAsFactors = FALSE)
genes <- read.csv(file.path(cellphone_data, "/gene_input.csv"), stringsAsFactors = FALSE)
proteins <- read.csv(file.path(cellphone_data, "/protein_input.csv"), stringsAsFactors = FALSE)

# subset the pbmc data to fewer cells to meet package requirements
RNA_features <- c(
    "CCL20", "CXCR3", "CCR6",
    "IL7", "IL7R", "IL2RG",
    "TGFB3", "TGFBR3",
    "ITGA6", "ITGB4", "NRG1"
)
TF_features <- c(
    "DBP",
    "FLI1", "ZNF431",
    "ZNF324", "CREM", "FOSL1"
)
name_features <- c(
    "CCL20", "CXCR3", "CCR6",
    "IL7", "IL7_receptor",
    "TGFB3", "TGFBR3",
    "integrin_a6b4_complex", "NRG1"
)

# From here, we start to diverge from the code for exported data,
# as some functions require multiple objects for testing
# Create 3 downsampled data sets with subsetted clusters and featuers
sub_celltypes <- c("CD8_T_cell", "CD14_monocyte", "B_cell")
n_cells <- 110
n_objects <- 3

set.seed(123)

# Get cell barcodes by cell type
cells_by_type <- setNames(lapply(sub_celltypes, function(clust) {
    colnames(pbmc)[pbmc$cell_type == clust]
}), sub_celltypes)

# Double check there are enough cells for non-overlapping samples
req_number <- n_cells * n_objects
too_small <- lengths(cells_by_type) < req_number
if (any(too_small)) {
    stop("Not enough cells in the following cell types for non-overlapping samples: ",
        toString(names(cells_by_type)[too_small]))
}

# For each cell type, sample and then divide into 3 non-overlapping groups:
subcells_by_type <- lapply(cells_by_type, sample, size = req_number, replace = FALSE)
cells_for_objects <- lapply(subcells_by_type, function(cells) {
    split(cells, ceiling(seq_along(cells) / n_cells))
})

# Double check for no overlap:
all_cells <- unlist(cells_for_objects, use.names = FALSE)
if (!identical(length(unique(all_cells)), length(all_cells))) {
    stop("There are overlapping cells across the objects.")
}

# Get cell barcodes for each object by combining across cell types
subset_cells <- lapply(seq_len(n_objects), function(i) {
    unlist(lapply(cells_for_objects, `[[`, i), use.names = FALSE)
})

# Extract counts and z-scores from SCE object
counts <- assay(pbmc, "counts")
z_scores <- t(scale(t(assay(pbmc, "logcounts"))))

# subset CellPhoneDB inputs
complexes_tiny <- complexes[complexes$complex_name %in% name_features, ]
genes_tiny <- genes[genes$gene_name %in% RNA_features, ]
proteins_tiny <- proteins[proteins$uniprot %in% genes_tiny$uniprot, ]
interactions_tiny <- interactions[
    (interactions$partner_a %in% proteins_tiny$uniprot | interactions$partner_a %in% complexes_tiny$complex_name) &
        (interactions$partner_b %in% proteins_tiny$uniprot | interactions$partner_b %in% complexes_tiny$complex_name),
]

# Format SCENIC inputs
auc <- t(auc)
rownames(auc) <- gsub("\\.\\.\\.$", "", rownames(auc))

regulons <- regulons[-1:-2, ]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
    "Annotation", "Context", "TargetGenes", "RankAtMax")
regulons_tiny <- regulons[regulons$TF %in% TF_features, ]

# Make rl_map
rl_map_tiny <- dominoSignal::create_rl_map_cellphonedb(
    genes = genes_tiny,
    proteins = proteins_tiny,
    interactions = interactions_tiny,
    complexes = complexes_tiny
)

# Get regulon list
regulon_list_tiny <- dominoSignal::create_regulon_list_scenic(
    regulons = regulons_tiny
)

# Get domino inputs (clusters, counts, z-scores)
get_inputs <- function(cells, n_cells, sub_celltypes, counts, z_scores, RNA_features) {
    clusts <- factor(rep(sub_celltypes, each = n_cells), levels = sub_celltypes)
    names(clusts) <- cells
    list(
        clusters = clusts,
        counts = counts[RNA_features, cells],
        z_scores = z_scores[RNA_features, cells],
        auc = auc[TF_features, cells]
    )
}

# Grab inputs for each object based on barcodes
dom_inputs <- lapply(subset_cells, get_inputs, n_cells = n_cells, sub_celltypes = sub_celltypes,
    counts = counts, z_scores = z_scores, RNA_features = RNA_features
)

# Create and build domino objects for each set of cells
pbmc_dom_tiny_list <- lapply(dom_inputs, function(inputs) {
    pbmc_dom_tiny <- dominoSignal::create_domino(
        rl_map = rl_map_tiny,
        features = inputs$auc,
        counts = inputs$counts,
        z_scores = inputs$z_scores,
        clusters = inputs$clusters,
        tf_targets = regulon_list_tiny,
        use_clusters = TRUE,
        use_complexes = TRUE,
        remove_rec_dropout = FALSE)
})

pbmc_dom_built_tiny_list <- lapply(pbmc_dom_tiny_list, dominoSignal::build_domino,
    min_tf_pval = 0.05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = 0.1,
    min_rec_percentage = 0.01
)

# Shorter names for development:
counts_tiny <- RNA_count_tiny
zscores_tiny <- RNA_zscore_tiny
tiny_created_dom1 <- pbmc_dom_tiny_list[[1]]
tiny_created_dom2 <- pbmc_dom_tiny_list[[2]]
tiny_created_dom3 <- pbmc_dom_tiny_list[[3]]
tiny_dom1 <- pbmc_dom_built_tiny_list[[1]]
tiny_dom2 <- pbmc_dom_built_tiny_list[[2]]
tiny_dom3 <- pbmc_dom_built_tiny_list[[3]]

# Save these to sysdata for lazy loading in development and testing; not exported for users
usethis::use_data(counts_tiny, zscores_tiny, clusters_tiny, rl_map_tiny, regulons_tiny,
    regulon_list_tiny, auc_tiny, complexes_tiny, genes_tiny, proteins_tiny, interactions_tiny,
    tiny_created_dom1, tiny_created_dom2, tiny_created_dom3, tiny_dom1, tiny_dom2, tiny_dom3,
    v0.2.1, internal = TRUE, overwrite = TRUE
)