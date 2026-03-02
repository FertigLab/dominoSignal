#' Convert genes using a table
#'
#' Takes a vector of gene inputs and a conversion table and returns a
#' converted gene table
#'
#' @param genes the genes to convert
#' @param from  gene symbol type of the input (ENSG, ENSMUSG, HGNC, MGI)
#' @param to    desired gene symbol type for the output (HGNC, MGI)
#' @param conversion_table a data frame with column names corresponding to gene symbol types (mm.ens, hs.ens, mgi, hgnc)
#' and rows corresponding to the gene symbols themselves
#' @return A data frame of genes with original and corresponding converted symbols
#' @keywords internal
#'
table_convert_genes <- function(genes, from, to, conversion_table) {
    # Check inputs:
    check_arg(genes, allow_class = c("character", "vector"))
    check_arg(from, allow_values = c("ENSMUSG", "ENSG", "MGI", "HGNC"))
    check_arg(to, allow_values = c("MGI", "HGNC"))
    check_arg(conversion_table, allow_class = "data.frame")
    stopifnot(`Conversion table must be provided with at least two of column names mm.ens, hs.ens, mgi and/or hgnc` =
            (length(which(colnames(conversion_table) %in% c(
                "mm.ens", "hs.ens", "mgi",
                "hgnc"
            ))) > 1))
    if (from == "ENSMUSG") {
        col1 <- conversion_table$mm.ens
    }
    if (from == "ENSG") {
        col1 <- conversion_table$hs.ens
    }
    if (from == "MGI") {
        col1 <- conversion_table$mgi
    }
    if (from == "HGNC") {
        col1 <- conversion_table$hgnc
    }
    if (to == "MGI") {
        col2 <- conversion_table$mgi
    }
    if (to == "HGNC") {
        col2 <- conversion_table$hgnc
    }
    genesV2 <- cbind(col1[which(col1 %in% genes)], col2[which(col1 %in% genes)])
    return(genesV2)
}

#' Use biomaRt to convert genes
#'
#' This function reads in a vector of genes and converts the genes to specified symbol type
#'
#' @param genes Vector of genes to convert.
#' @param from Format of gene input (ENSMUSG, ENSG, MGI, or HGNC)
#' @param to Format of gene output (MGI or HGNC)
#' @param host Host to connect to. Defaults to https://www.ensembl.org following the useMart default,
#'   but can be changed to archived hosts if useMart fails to connect.
#' @importFrom biomaRt useMart getLDS
#' @return A data frame with input genes as column 1 and converted genes as column 2
#' @keywords internal
#'
convert_genes <- function(
    genes, from = c("ENSMUSG", "ENSG", "MGI", "HGNC"), to = c("MGI", "HGNC"),
    host = "https://www.ensembl.org"
) {
    # Check inputs:
    check_arg(genes, allow_class = c("character", "vector"))
    check_arg(from, allow_values = c("ENSMUSG", "ENSG", "MGI", "HGNC"))
    check_arg(to, allow_values = c("MGI", "HGNC"))
    check_arg(host, allow_class = "character", allow_len = 1)

    if (from == "ENSMUSG") {
        srcMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
        sourceAtts <- "ensembl_gene_id"
    }
    if (from == "ENSG") {
        srcMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
        sourceAtts <- "ensembl_gene_id"
    }
    if (from == "MGI") {
        srcMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
        sourceAtts <- "mgi_symbol"
    }
    if (from == "HGNC") {
        srcMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
        sourceAtts <- "hgnc_symbol"
    }
    if (to == "MGI") {
        tarMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
        tarAtts <- "mgi_symbol"
    }
    if (to == "HGNC") {
        tarMart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
        tarAtts <- "hgnc_symbol"
    }
    genesV2 <- getLDS(
        attributes = sourceAtts, filters = sourceAtts, values = genes, mart = srcMart,
        attributesL = tarAtts, martL = tarMart, uniqueRows = FALSE
    )
    return(genesV2)
}
