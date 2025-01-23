

#' Generate Expression Bootstrap
#' 
#' Sample cell expression profiles, including transcription factor scores and cluster assignments, with replacement to generate a bootstrap for use in comparison of differential signaling between datasets without sample replicate annotation. Bootstraps will include as many cells from each cluster as the parent data set used to generate the bootstrap.
#' 
#' @param counts A matrix of RNA expression counts
#' @param normalized_expression A matrix of normalized RNA expression counts
#' @param clusters A factor of cell assignments to cell type clusters
#' @param TF_scores A matrix of scores representing TF activity in cells
#' @return A list containing expression values for the derived bootstrap:
#' \itemize{
#'  \item{'counts'} : A matrix of RNA expression counts
#'  \item{'scaled_expression'} : A matrix of scaled normalized RNA expression counts centered on 0
#'  \item{'clusters'} : A factor of cell assignments to cell type clusters
#'  \item{'TF_scores'} : A matrix of scores representing TF activity in cells
#' }
#' @export
#' 

generate_expression_bootstrap <- function(counts, normalized_expression, clusters, TF_scores){
  if(!is(clusters, "factor")){
    cluster_lvls <- levels(as.factor(clusters))
  } else {
    cluster_lvls <- levels(clusters)
  }
  
  # select cells for the bootstrap by sampling sample names with replacement
  cell_ids <- c()
  for(j in seq_along(cluster_lvls)){
    cl <- cluster_lvls[j]
    cl_bool <- clusters == cl
    cl_count <- sum(cl_bool)
    cl_ids <- colnames(norm_expr)[cl_bool]
    
    boot_ids <- sample(cl_ids, cl_count, replace = TRUE)
    cell_ids <- c(cell_ids, boot_ids)
  }
  
  # create bootstrap matrices
  counts_i <- counts[, cell_ids]
  logcounts <- norm_expr[rowSums(norm_expr) != 0, cell_ids]
  if(!is(logcounts, "matrix")){
    logcounts <- as.matrix(logcounts)
  }
  z_scores_i <- t(scale(t(logcounts)))
  clusters_i <- clusters[cell_ids] %>% factor(., levels = cluster_lvls)
  TF_scores_i <- TF_scores[,cell_ids]
  
  res <- list(
    "counts" = counts_i,
    "scaled_expression" = z_scores_i,
    "clusters" = clusters_i,
    "TF_scores" = TF_scores_i
  )
  return(res)
}

