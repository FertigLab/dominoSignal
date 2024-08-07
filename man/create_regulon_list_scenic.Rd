% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_fxns.R
\name{create_regulon_list_scenic}
\alias{create_regulon_list_scenic}
\title{Create a list of genes in regulons inferred by SCENIC}
\usage{
create_regulon_list_scenic(regulons)
}
\arguments{
\item{regulons}{Data frame or file path to the table of the output of the ctx function from pySCENIC}
}
\value{
A list where names are transcription factors and the stored values are character vectors of genes in the inferred regulons
}
\description{
Generates a list of transcription factors and the genes targeted by the transcription factor as part of their regulon inferred by pySCENIC
}
\examples{
data(SCENIC)
regulon_list_tiny <- create_regulon_list_scenic(regulons = SCENIC$regulons_tiny)

}
