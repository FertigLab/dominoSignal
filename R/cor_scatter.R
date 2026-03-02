#' Create a correlation plot between TF and receptor
#'
#' Create a correlation plot between transcription factor activation score and receptor expression
#'
#' @param dom Domino object with network built ([build_domino()])
#' @param tf Target TF for plottting AUC score
#' @param rec Target receptor for plotting expression
#' @param remove_rec_dropout Whether to remove cells with zero expression for plot.
#'  This should match the same setting as in [build_domino()].
#' @param ... Other parameters to pass to [ggpubr::ggscatter()].
#' @return A ggplot scatter plot rendered in the active graphics device
#' @export cor_scatter
#' @examples
#' example(build_domino, echo = FALSE)
#' cor_scatter(pbmc_dom_built_tiny, "FLI1","CXCR3")
#'
cor_scatter <- function(dom, tf, rec, remove_rec_dropout = TRUE, ...) {
    if (remove_rec_dropout) {
        keep_id <- which(dom@counts[rec, ] > 0)
        rec_z_scores <- dom@z_scores[rec, keep_id]
        tar_tf_scores <- dom@features[tf, keep_id]
    } else {
        rec_z_scores <- dom@z_scores[rec, ]
        tar_tf_scores <- dom@features[tf, ]
    }
    dat <- data.frame(rec = rec_z_scores, tf = tar_tf_scores)
    ggscatter(dat, x = "rec", y = "tf", add = "reg.line", conf.int = FALSE, cor.coef = FALSE,
        cor.method = "pearson", xlab = rec, ylab = tf, size = 0.25, ...)
}
