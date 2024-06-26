% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dom_counts}
\alias{dom_counts}
\title{Access counts}
\usage{
dom_counts(dom)
}
\arguments{
\item{dom}{A domino object that has been created with \code{\link[=create_domino]{create_domino()}}}
}
\value{
A matrix containing the gene expression values for each gene (row) by cell (column)
}
\description{
A function to pull gene expression from a domino object
}
\examples{
counts <- dom_counts(dominoSignal:::pbmc_dom_built_tiny)

}
