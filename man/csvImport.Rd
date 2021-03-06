% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csvImport.R
\name{csvImport}
\alias{csvImport}
\title{support for importing separate CSV with fat format: one record per individual}
\usage{
csvImport(gpath, ipath, spath, times = c(0, 10, 20, 30, 60, 90, 120),
  gfilt = force, ifilt = force, sfilt = force)
}
\arguments{
\item{gpath}{character pathname of glucose data}

\item{ipath}{character pathname of insulin data}

\item{spath}{character pathname of sample-level data}

\item{times}{numeric vector of observation times in minutes}

\item{gfilt}{function that operates on glucose data accepting and returning data.frame as returned by read.csv}

\item{ifilt}{function that operates on insulin data accepting and returning data.frame as returned by read.csv}

\item{sfilt}{function that operates on sample data accepting and returning data.frame as returned by read.csv}
}
\description{
support for importing separate CSV with fat format: one record per individual
}
\examples{
pref = system.file("csv_example", package="ogttMetrics")
gpath = paste0(pref, "/glucBase.csv")
ipath = paste0(pref, "/insBase.csv")
spath = paste0(pref, "/sampBase.csv")
democ = csvImport(gpath, ipath, spath)
democ
}
