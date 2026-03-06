# generate tiny objects for examples in functions and for users

library(SingleCellExperiment)
library(dominoSignal)

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

# Create 3 downsampled data sets with subsetted clusters and featuers
sub_celltypes <- c("CD8_T_cell", "CD14_monocyte", "B_cell")
n_cells <- 120

set.seed(123)
for (i in seq_along(cell_types_dwn)) {
    cell <- cell_types_dwn[i]
    cell_barcodes <- colnames(pbmc)[pbmc$cell_type == cell]
    dwn_barcodes <- sample(cell_barcodes, n, replace = FALSE)
    cell_list[[cell]] <- dwn_barcodes
}
barcodes_dwn <- unlist(cell_list)
clusters_tiny <- factor(rep(names(cell_list), lengths(cell_list)))
names(clusters_tiny) <- barcodes_dwn
# Extract counts and z-scores from SCE object

counts <- assay(pbmc, "counts")
z_scores <- t(scale(t(assay(pbmc, "logcounts"))))

RNA_count_tiny <- counts[
    rownames(assay(pbmc, "counts")) %in% RNA_features,
    colnames(assay(pbmc, "counts")) %in% barcodes_dwn
]
RNA_zscore_tiny <- z_scores[
    rownames(assay(pbmc, "logcounts")) %in% RNA_features,
    colnames(assay(pbmc, "logcounts")) %in% barcodes_dwn
]


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
auc_tiny <- auc[TF_features, barcodes_dwn]

regulons <- regulons[-1:-2, ]
colnames(regulons) <- c(
    "TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
    "Annotation", "Context", "TargetGenes", "RankAtMax"
)
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

pbmc_dom_tiny <- dominoSignal::create_domino(
    rl_map = rl_map_tiny,
    features = auc_tiny,
    counts = RNA_count_tiny,
    z_scores = RNA_zscore_tiny,
    clusters = clusters_tiny,
    tf_targets = regulon_list_tiny,
    use_clusters = TRUE,
    use_complexes = TRUE,
    remove_rec_dropout = FALSE
)

pbmc_dom_built_tiny <- dominoSignal::build_domino(
    dom = pbmc_dom_tiny,
    min_tf_pval = 0.05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = 0.1,
    min_rec_percentage = 0.01
)

# Make a mock linkage summary as well:
linkage_sum_tiny <- new("linkage_summary",
    subject_meta = data.frame(
        "subject_names" = paste0("P", 1:6),
        "group" = c(rep("G1", 3), rep("G2", 3))
    ),
    subject_names = factor(
        paste0("P", 1:6),
        levels = paste0("P", 1:6)
    ),
    subject_linkages = list(
        "P1" = list(
            "C1" = list(
                "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                "rec" = c("R1", "R2", "R3", "R4"),
                "incoming_lig" = c("L1", "L2", "L3", "L4"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
            ),
            "C2" = list(
                "tfs" = c("TF2", "TF3", "TF4"),
                "rec" = c("R2", "R3", "R4"),
                "incoming_lig" = c("L2", "L3", "L4"),
                "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
            )
        ),
        "P2" = list(
            "C1" = list(
                "tfs" = c("TF1", "TF2"),
                "rec" = c("R1", "R2"),
                "incoming_lig" = c("L1", "L2"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2")
            ),
            "C2" = list(
                "tfs" = c("TF3", "TF4"),
                "rec" = c("R3", "R4"),
                "incoming_lig" = c("L3", "L4"),
                "tfs_rec" = c("TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R3 <- L3", "R4 <- L4")
            )
        ),
        "P3" = list(
            "C1" = list(
                "tfs" = c("TF1", "TF2"),
                "rec" = c("R1", "R2"),
                "incoming_lig" = c("L1", "L2"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2")
            ),
            "C2" = list(
                "tfs" = "TF3",
                "rec" = "R3",
                "incoming_lig" = "L3",
                "tfs_rec" = "TF3 <- R3",
                "rec_lig" = "R3 <- L3"
            )
        ),
        "P4" = list(
            "C1" = list(
                "tfs" = c("TF2", "TF3", "TF4"),
                "rec" = c("R2", "R3", "R4"),
                "incoming_lig" = c("L2", "L3", "L4"),
                "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
            ),
            "C2" = list(
                "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                "rec" = c("R1", "R2", "R3", "R4"),
                "incoming_lig" = c("L1", "L2", "L3", "L4"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
            )
        ),
        "P5" = list(
            "C1" = list(
                "tfs" = "TF3",
                "rec" = "R3",
                "incoming_lig" = "L3",
                "tfs_rec" = "TF3 <- R3",
                "rec_lig" = "R3 <- L3"
            ),
            "C2" = list(
                "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                "rec" = c("R1", "R2", "R3", "R4"),
                "incoming_lig" = c("L1", "L2", "L3", "L4"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
            )
        ),
        "P6" = list(
            "C1" = list(
                "tfs" = NULL,
                "rec" = NULL,
                "incoming_lig" = NULL,
                "tfs_rec" = NULL,
                "rec_lig" = NULL
            ),
            "C2" = list(
                "tfs" = c("TF1", "TF2", "TF3"),
                "rec" = c("R1", "R2", "R3"),
                "incoming_lig" = c("L1", "L2", "L3"),
                "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3"),
                "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3")
            )
        )
    )
)

# save all data to be used in examples and for users to explore (exported)
CellPhoneDB <- list(
    complexes_tiny = complexes_tiny,
    genes_tiny = genes_tiny,
    proteins_tiny = proteins_tiny,
    interactions_tiny = interactions_tiny,
    rl_map_tiny = rl_map_tiny
)
save(CellPhoneDB, file = file.path("data", "CellPhoneDB.RData"))

SCENIC <- list(
    auc_tiny = auc_tiny,
    regulons_tiny = regulons_tiny,
    regulon_list_tiny = regulon_list_tiny
)
save(SCENIC, file = file.path("data", "SCENIC.RData"))

PBMC <- list(
    count_tiny = RNA_count_tiny,
    zscore_tiny = RNA_zscore_tiny,
    clusters_tiny = clusters_tiny
)
save(PBMC, file = file.path("data", "PBMC.RData"))

DominoObjects <- list(
    dom_tiny = pbmc_dom_tiny,
    built_dom_tiny = pbmc_dom_built_tiny
)
save(DominoObjects, file = file.path("data", "DominoObjects.RData"))

LinkageSummary <- linkage_sum_tiny
save(LinkageSummary, file = file.path("data", "LinkageSummary.RData"))
