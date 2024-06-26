% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_fxns.R
\name{lc}
\alias{lc}
\title{Pulls all items from a list pooled into a single vector}
\usage{
lc(list, list_names)
}
\arguments{
\item{list}{List to pull items from}

\item{list_names}{Names of items in list to pool}
}
\value{
A vector contaning all items in the list by list_names
}
\description{
Helper function to convert from a nested series of lists to a single vector.
}
\keyword{internal}
