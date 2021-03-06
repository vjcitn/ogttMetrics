% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIexplorer.R
\name{SIexplorer}
\alias{SIexplorer}
\title{interactive multipanel visualization for ogttCohort instance}
\usage{
SIexplorer(oc = obaSamp, winsorizeSI = TRUE, ...)
}
\arguments{
\item{oc}{ogttCohort instance}

\item{winsorizeSI}{if TRUE, move negative estimates of SI to smallest positive value}

\item{\dots}{passed to \code{\link{minmodByID}}}
}
\description{
interactive multipanel visualization for ogttCohort instance
}
\examples{
if (interactive()) {
  if (options()$example.ask) stop("must set options(example.ask=FALSE) before running example")
  data(obaSamp)
  SIexplorer(obaSamp)
}
}
