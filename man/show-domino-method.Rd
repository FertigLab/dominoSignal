% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_definitions.R
\name{show,domino-method}
\alias{show,domino-method}
\title{Show domino object information}
\usage{
\S4method{show}{domino}(object)
}
\arguments{
\item{object}{A domino object}
}
\value{
A printed description of cell numbers and clusters in the object
}
\description{
Shows content overview of domino object
}
\examples{
example(build_domino, echo = FALSE)
show(pbmc_dom_built_tiny)

}
