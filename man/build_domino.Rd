% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_fxns.R
\name{build_domino}
\alias{build_domino}
\title{Calculate a signaling network for a domino object}
\usage{
build_domino(
  dom,
  max_tf_per_clust = 5,
  min_tf_pval = 0.01,
  max_rec_per_tf = 5,
  rec_tf_cor_threshold = 0.15,
  min_rec_percentage = 0.1
)
}
\arguments{
\item{dom}{Domino object from create_domino.}

\item{max_tf_per_clust}{Maximum number of transcription factors called active in a cluster.}

\item{min_tf_pval}{Minimum p-value from differential feature score test to call a transcription factor active in a cluster.}

\item{max_rec_per_tf}{Maximum number of receptors to link to each transcription factor.}

\item{rec_tf_cor_threshold}{Minimum pearson correlation used to consider a receptor linked with a transcription factor. Increasing this will decrease the number of receptors linked to each transcription factor.}

\item{min_rec_percentage}{Minimum percentage of cells in cluster expressing a receptor for the receptor to be linked to trancription factors in that cluster.}
}
\value{
A domino object with a signaling network built
}
\description{
This function calculates a signaling network. It requires a domino object
preprocessed from create_domino and returns a domino object prepared for
plotting with the various plotting functions in this package.
}
\examples{
pbmc_dom_tiny_built <- build_domino(
 dom = dominoSignal:::pbmc_dom_tiny, min_tf_pval = .001, max_tf_per_clust = 25,
 max_rec_per_tf = 25, rec_tf_cor_threshold = .25, min_rec_percentage = 0.1
)

}
