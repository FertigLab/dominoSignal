#' @import biomaRt
#' @import stats
#' @importFrom utils read.csv
#' @importFrom Matrix rowSums
#' @import methods
#'
NULL

#' Create a receptor - ligand map from a CellPhoneDB signaling database
#'
#' Generates a data frame of ligand-receptor interactions from a CellPhoneDB database annotating the genes
#'   encoding the interacting ligands and receptors to be queried in transcriptomic data.
#'
#' @param genes data frame or file path to table of gene names in uniprot, hgnc_symbol, or ensembl format in
#'   CellPhoneDB database format
#' @param proteins data frame or file path to table of protein features in CellPhoneDB format
#' @param interactions data frame or file path to table of protein-protein interactions in CellPhoneDB format
#' @param complexes optional: data frame or file path to table of protein complexes in CellPhoneDB format
#' @param database_name name of the database being used, stored in output
#' @param gene_conv a tuple of (from, to) or (source, target) if gene conversion to orthologs is desired;
#'   options are ENSMUSG, ENSG, MGI, or HGNC
#' @param gene_conv_host host for conversion; default ensembl, could also use mirrors if desired
#' @param alternate_convert boolean if you would like to use a non-ensembl method of conversion
#'   (must supply table; not recommended, use only if ensembl is down)
#' @param alternate_convert_table a data frame with column names corresponding to gene symbol types
#'   (mm.ens, hs.ens, mgi, hgnc) and rows corresponding to the gene symbols themselves for use with
#'   alternate_convert = TRUE
#' @return Data frame where each row describes a possible receptor-ligand interaction
#' @export create_rl_map_cellphonedb
#' @examples
#' data(CellPhoneDB)
#' rl_map_tiny <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
#'  proteins = CellPhoneDB$proteins_tiny,
#'  interactions = CellPhoneDB$interactions_tiny,
#'  complexes = CellPhoneDB$complexes_tiny)
#' 
#' \dontrun{
#' rl_map_tiny_conv <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
#'   proteins = CellPhoneDB$proteins_tiny,
#'   interactions = CellPhoneDB$interactions_tiny,
#'   gene_conv = c("HGNC", "MGI"),
#'   gene_conv_host = "https://beta.ensembl.org")
#' }
#' 
#' # Using alternate conversion table instead of biomaRt
#' ortho_table <- data.frame(
#'   hs.ens = c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804"),
#'   hgnc = c("MT-ND1", "MT-ND2", "MT-CO1"),
#'   mm.ens = c("ENSMUSG00000064341", "ENSMUSG00000064345", "ENSMUSG00000064351"),
#'   mgi = c("mt-Nd1", "mt-Nd2", "mt-Co1"))
#' 
#' rl_map_tiny_alt <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
#'   proteins = CellPhoneDB$proteins_tiny,
#'   interactions = CellPhoneDB$interactions_tiny,
#'   complexes = CellPhoneDB$complexes_tiny, gene_conv = c("ENSG", "MGI"),
#'   alternate_convert = TRUE, alternate_convert_table = ortho_table)
#' 
create_rl_map_cellphonedb <- function(
    genes, proteins, interactions, complexes = NULL, database_name = "CellPhoneDB",
    gene_conv = NULL, gene_conv_host = "https://www.ensembl.org", alternate_convert = FALSE,
    alternate_convert_table = NULL
) {
    # Check input structures:
    check_arg(genes, c("character", "data.frame"))
    check_arg(proteins, c("character", "data.frame"))
    check_arg(interactions, c("character", "data.frame"))
    check_arg(complexes, c("character", "data.frame", "NULL"))
    check_arg(database_name, "character", allow_len = 1)
    check_arg(gene_conv, c("NULL", "character"), allow_len = c(0, 2))
    check_arg(gene_conv_host, "character", allow_len = 1)

    # Read in files if needed:
    genes <- read_if_char(genes)
    proteins <- read_if_char(proteins)
    interactions <- read_if_char(interactions)
    complexes <- read_if_char(complexes)

    # replace empty cells in columns annotating gene properties with 'False' There are some
    # unannotated genes in database v2.0 that seem to have been fixed in v4.0
    gene_features <- c(
        "transmembrane", "peripheral", "secreted", "secreted_highlight", "receptor",
        "integrin", "other"
    )
    proteins[!nzchar(proteins$receptor, keepNA = TRUE), colnames(proteins) %in% gene_features] <- "False"

    # change cases of True/False syntax from Python to TRUE/FALSE R syntax
    genes <- conv_py_bools(genes)
    proteins <- conv_py_bools(proteins)
    interactions <- conv_py_bools(interactions)
    complexes <- conv_py_bools(complexes)

    # gene conversions
    if (!is.null(gene_conv) && !identical(gene_conv[1], gene_conv[2])) {
        # obtain conversion dictionary
        if (alternate_convert) {
            conv_dict <- table_convert_genes(genes$gene_name,
                from = gene_conv[1], to = gene_conv[2],
                alternate_convert_table
            )
        } else {
            conv_dict <- convert_genes(genes$gene_name, from = gene_conv[1], to = gene_conv[2], host = gene_conv_host)
        }
        # column 1 is the source gene names used by the reference data base column 2 is the
        # orthologous gene names for the organism to which the reference is being converted
    }
    # Step through the interactions and build rl connections.
    rl_map <- NULL
    for (i in seq_len(nrow(interactions))) {
        inter <- interactions[i, ]
        partner_a <- inter[["partner_a"]]
        partner_b <- inter[["partner_b"]]
        conversion_flag <- list()
        # features of partner_a
        a_features <- list()
        if (partner_a %in% complexes[["complex_name"]]) {
            complex_a <- complexes[complexes[["complex_name"]] == partner_a, ]
            component_a <- as.character(complex_a[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")])
            component_a <- component_a[nzchar(component_a, keepNA = TRUE)]
            a_features[["uniprot_A"]] <- paste(component_a, collapse = ",")
            gene_a <- vapply(component_a, FUN.VALUE = character(1), FUN = function(x) {
                g <- unique(genes[genes[["uniprot"]] == x, "gene_name"])
                if (!is.null(gene_conv) && !identical(gene_conv[1], gene_conv[2])) {
                    # if the original gene trying to be converted is not in the gene dictionary the
                    # interaction is not included in the final rl_map
                    if (sum(g %in% conv_dict[, 1]) < length(g)) {
                        for (gn in g) {
                            conversion_flag[[gn]] <- TRUE
                        }
                    } else {
                        g <- paste(unique(conv_dict[conv_dict[, 1] %in% g, 2]), collapse = ";")
                    }
                }
                # if multiple genes are annotated for the uniprot ID, use only the first unique instance
                if (length(g) == 1) {
                    res <- g
                } else {
                    res <- g[1]
                    g_col <- toString(g)
                    message(
                        component_a, " has multiple encoding gene mapped in genes table.\n",
                        g_col, "\n",
                        "The first mapping gene is used: ", res
                    )
                }
                return(res)
            })
            a_features[["gene_A"]] <- paste(gene_a, collapse = ",")
            # annotation as a receptor or ligand is based on the annotation of the complex
            a_features[["type_A"]] <- ifelse(complex_a[["receptor"]], "R", "L")
            # replace any spaces in the partner name with an underscore
            a_features[["name_A"]] <- gsub(" ", "_", partner_a, fixed = TRUE)
        } else if (partner_a %in% proteins[["uniprot"]]) {
            protein_a <- proteins[proteins[["uniprot"]] == partner_a, ]
            component_a <- protein_a[["uniprot"]]
            a_features[["uniprot_A"]] <- component_a
            gene_a <- unique(genes[genes[["uniprot"]] == component_a, "gene_name"])
            if (!is.null(gene_conv) && !identical(gene_conv[1], gene_conv[2])) {
                # if the original gene trying to be converted is not in the gene dictionary the
                # interaction is not included in the final rl_map
                if (sum(gene_a %in% conv_dict[, 1]) < length(gene_a)) {
                    for (gn in gene_a) {
                        conversion_flag[[gn]] <- TRUE
                    }
                } else {
                    gene_a <- unique(conv_dict[conv_dict[, 1] %in% gene_a, 2])
                }
            }
            gene_a <- paste(gene_a, collapse = ";")
            a_features[["gene_A"]] <- gene_a
            a_features[["type_A"]] <- ifelse(protein_a[["receptor"]], "R", "L")
            a_features[["name_A"]] <- gene_a
        } else {
            next
        }
        if (length(conversion_flag)) {
            message(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
            message(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
            next
        }
        a_df <- as.data.frame(a_features)
        # features of partner_b
        b_features <- list()
        if (partner_b %in% complexes[["complex_name"]]) {
            complex_b <- complexes[complexes[["complex_name"]] == partner_b, ]
            component_b <- as.character(complex_b[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")])
            component_b <- component_b[nzchar(component_b, keepNA = TRUE)]
            b_features[["uniprot_B"]] <- paste(component_b, collapse = ",")
            gene_b <- vapply(component_b, FUN.VALUE = character(1), FUN = function(x) {
                g <- unique(genes[genes[["uniprot"]] == x, "gene_name"])
                if (!is.null(gene_conv) && !identical(gene_conv[1], gene_conv[2])) {
                    # if the original gene trying to be converted is not in the gene dictionary the
                    # interaction is not included in the final rl_map
                    if (sum(g %in% conv_dict[, 1]) < length(g)) {
                        for (gn in g) {
                            conversion_flag[[gn]] <- TRUE
                        }
                    } else {
                        g <- paste(unique(conv_dict[conv_dict[, 1] %in% g, 2]), collapse = ";")
                    }
                }
                # if multiple genes are annotated for the uniprot ID, use only the first unique instance
                if (length(g) == 1) {
                    res <- g
                } else {
                    res <- g[1]
                    g_col <- toString(g)
                    message(
                        component_a, " has multiple encoding gene mapped in genes table.\n",
                        g_col, "\n",
                        "The first mapping gene is used: ", res
                    )
                }
                return(res)
            })
            b_features[["gene_B"]] <- paste(gene_b, collapse = ",")
            # annotation as a receptor or ligand is based on the annotation of the complex
            b_features[["type_B"]] <- ifelse(complex_b[["receptor"]], "R", "L")
            # replace any spaces in the partner name with an underscore
            b_features[["name_B"]] <- gsub(" ", "_", partner_b, fixed = TRUE)
        } else if (partner_b %in% proteins[["uniprot"]]) {
            protein_b <- proteins[proteins[["uniprot"]] == partner_b, ]
            component_b <- protein_b[["uniprot"]]
            b_features[["uniprot_B"]] <- component_b
            gene_b <- unique(genes[genes[["uniprot"]] == component_b, "gene_name"])
            if (!is.null(gene_conv) && !identical(gene_conv[1], gene_conv[2])) {
                # if the original gene trying to be converted is not in the gene dictionary the
                # interaction is not included in the final rl_map
                if (sum(gene_b %in% conv_dict[, 1]) < length(gene_b)) {
                    for (gn in gene_b) {
                        conversion_flag[[gn]] <- TRUE
                    }
                } else {
                    gene_b <- unique(conv_dict[conv_dict[, 1] %in% gene_b, 2])
                }
            }
            gene_a <- paste(gene_b, collapse = ";")
            b_features[["gene_B"]] <- gene_b
            b_features[["type_B"]] <- ifelse(protein_b[["receptor"]], "R", "L")
            b_features[["name_B"]] <- gsub(" ", "_", gene_b, fixed = TRUE)
        } else {
            next
        }
        if (length(conversion_flag)) {
            message(paste("No gene orthologs found for:", names(conversion_flag), collapse = " "))
            message(paste("Skipping interaction:", partner_a, partner_b, collapse = " "))
            next
        }
        b_df <- as.data.frame(b_features)
        i_features <- cbind(a_df, b_df)
        i_features[["int_pair"]] <- paste(i_features[["name_A"]], i_features[["name_B"]], sep = " & ")
        i_features[["annotation_strategy"]] <- inter[["annotation_strategy"]]
        i_features[["source"]] <- inter[["source"]]
        i_features[["database_name"]] <- database_name
        rl_map <- rbind(i_features, rl_map)
    }
    # exclude rows without receptor-ligand interactions
    rl_map <- rl_map[!(rl_map$type_A == "R" & rl_map$type_B == "R") & !(rl_map$type_A == "L" & rl_map$type_B == "L"), ]
    # specify column order
    rl_map <- rl_map[, c(
        "int_pair", "name_A", "uniprot_A", "gene_A", "type_A", "name_B", "uniprot_B",
        "gene_B", "type_B", "annotation_strategy", "source", "database_name"
    )]
    return(rl_map)
}