% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_arg}
\alias{check_arg}
\title{Check input arguments}
\usage{
check_arg(
  arg,
  allow_class = NULL,
  allow_len = NULL,
  allow_range = NULL,
  allow_values = NULL,
  need_vars = c(NULL),
  need_colnames = FALSE,
  need_rownames = FALSE,
  need_names = FALSE
)
}
\arguments{
\item{arg}{the argument to check}

\item{allow_class}{vector of allowed classes}

\item{allow_len}{vector of allowed lengths}

\item{allow_range}{range of minimum and maximum values i.e. c(1, 5)}

\item{allow_values}{vector of allowed values}

\item{need_vars}{vector of required variables}

\item{need_colnames}{vogical for whether colnames are required}

\item{need_rownames}{logical for whether rownames are required}

\item{need_names}{logical for whether names are required}
}
\value{
Logical indicating whether the argument meets the requirements
}
\description{
Accepts an object and rules to check against; stops if requirements are not met
}
\keyword{internal}
