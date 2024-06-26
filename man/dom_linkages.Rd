% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dom_linkages}
\alias{dom_linkages}
\title{Access linkages}
\usage{
dom_linkages(
  dom,
  link_type = c("complexes", "receptor-ligand", "tf-target", "tf-receptor", "receptor",
    "incoming-ligand"),
  by_cluster = FALSE
)
}
\arguments{
\item{dom}{A domino object that has been created with \code{\link[=create_domino]{create_domino()}}}

\item{link_type}{One value (out of "complexes", "receptor-ligand",
"tf-target", "tf-receptor", "receptor", "incoming-ligand") used
to select the desired type of linkage}

\item{by_cluster}{A boolean to indicate whether the linkages should be returned overall or by cluster}
}
\value{
A list containing linkages between some combination of receptors, ligands, transcription factors, and clusters
}
\description{
A function to pull linkages from a domino object
}
\examples{
complexes <- dom_linkages(dominoSignal:::pbmc_dom_built_tiny, "complexes")
tf_rec_by_cluster <- dom_linkages(dominoSignal:::pbmc_dom_built_tiny, "tf-receptor", TRUE)

}
