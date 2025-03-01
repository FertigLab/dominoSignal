% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fxns.R
\name{obtain_circos_expression}
\alias{obtain_circos_expression}
\title{Obtain Circos Expression}
\usage{
obtain_circos_expression(
  dom,
  receptor,
  ligands,
  ligand_expression_threshold = 0.01,
  cell_idents = NULL
)
}
\arguments{
\item{dom}{Domino object that has undergone network building with build_domino()}

\item{receptor}{Name of a receptor active in at least one cell type in the domino object}

\item{ligands}{Character vector of ligands capable of interaction with the receptor}

\item{ligand_expression_threshold}{Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor}

\item{cell_idents}{Vector of cell types from cluster assignments in the domino object to be included in the plot.}
}
\value{
a data frame where each row describes plotting parameters of ligand-receptor interactions to pass to render_circos_ligand_receptor()
}
\description{
Pull expression data from a domino object and format for plotting as a receptor-oriented circos plot.
}
\examples{
example(build_domino, echo = FALSE)
#basic usage
obtain_circos_expression(pbmc_dom_built_tiny, receptor = "CXCR3", ligands = "CCL20")

}
