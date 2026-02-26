#' Create a list of genes in regulons inferred by SCENIC
#'
#' Generates a list of transcription factors and the genes targeted by the transcription factor as part of their regulon inferred by pySCENIC
#'
#' @param regulons Data frame or file path to the table of the output of the ctx function from pySCENIC
#' @return A list where names are transcription factors and the stored values are character vectors of genes in the inferred regulons
#' @export create_regulon_list_scenic
#' @examples
#' data(SCENIC)
#' regulon_list_tiny <- create_regulon_list_scenic(regulons = SCENIC$regulons_tiny)
#'
create_regulon_list_scenic <- function(regulons) {
  if (is(regulons, "character")) {
    regulons <- read.csv(regulons)
  }
  TFS <- unique(regulons[["TF"]])
  TF_targets <- lapply(TFS, function(tf) {
    regulons_small <- regulons[regulons[["TF"]] == tf, ]
    targets <- regulons_small[["TargetGenes"]]
    target_genes <- lapply(targets, function(x) {
      split_targs <- unlist(strsplit(x, ""))
      split_targs <- split_targs[seq(2, length(split_targs))]
      split_targs <- split_targs[seq(1, length(split_targs) - 1)]
      split_targs <- paste(split_targs, collapse = "")
      split_targs <- unlist(strsplit(split_targs, "['), (']"))
      split_targs_pre <- split_targs[split_targs != ""]
      split_targs_post <- split_targs_pre[seq(1, length(split_targs_pre), 2)]
      return(split_targs_post)
    })
    return(unique(unlist(target_genes)))
  })
  names(TF_targets) <- TFS
  return(TF_targets)
}
