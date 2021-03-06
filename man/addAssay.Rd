% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allS4.R
\name{addAssay}
\alias{addAssay}
\title{add a matrix-valued assay to a MultiAssayExperiment instance}
\usage{
addAssay(mae, mat, assayname, allowGaps = FALSE)
}
\arguments{
\item{mae}{MultiAssayExperiment instance}

\item{mat}{a matrix with rownames assigned; if colnames absent will be taken from rownames(colData(mae))}

\item{assayname}{character string used to name assay in ExperimentList of mae}

\item{allowGaps}{logical, should addAssay fail if there are gaps between sample labels for new assay and rownames of colData}
}
\description{
add a matrix-valued assay to a MultiAssayExperiment instance
}
\examples{
data(obaSamp)
lit = subsetByAssay(obaSamp, c("glucose", "insulin"))
d = dim(ee <- experiments(lit)$glucose)
dum = matrix(rnorm(prod(d)), nrow=d[1])
rownames(dum) = rownames(ee)
colnames(dum) = colnames(ee)
addAssay(lit, dum, "dum")
}
