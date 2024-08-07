% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_definitions.R
\docType{class}
\name{domino-class}
\alias{domino-class}
\alias{domino}
\title{The domino class}
\value{
An instance of class \code{domino }
}
\description{
The domino class contains all information necessary to calculate receptor-ligand
signaling. It contains z-scored expression, cell cluster labels, feature values,
and a referenced receptor-ligand database formatted as a receptor-ligand map.
Calculated intermediate values are also stored.
}
\section{Slots}{

\describe{
\item{\code{db_info}}{List of data sets from ligand - receptor database}

\item{\code{counts}}{Raw count gene expression data}

\item{\code{z_scores}}{Matrix of z-scored expression data with cells as columns}

\item{\code{clusters}}{Named factor with cluster identity of each cell}

\item{\code{features}}{Matrix of features (TFs) to correlate receptor - ligand expression with. Cells are columns and features are rows.}

\item{\code{cor}}{Correlation matrix of receptor expression to features.}

\item{\code{linkages}}{List of lists containing info linking cluster->tf->rec->lig}

\item{\code{clust_de}}{Data frame containing differential expression results for features by cluster.}

\item{\code{misc}}{List of miscellaneous info pertaining to run parameters etc.}

\item{\code{cl_signaling_matrices}}{Incoming signaling matrix for each cluster}

\item{\code{signaling}}{Signaling matrix between all clusters.}
}}

