
matsuda120 = function( gluc, ins, times ) {
#
# http://mmatsuda.diabetes-smc.jp/MIndex.html has link to excel macro
#
 stopifnot(length(gluc)==length(ins), 
   length(gluc)==length(times))
 timesUsed = c(0, 30, 60, 90, 120)
 keep = which(times %in% timesUsed)
 stopifnot(all.equal(times[keep], timesUsed))
 gluc = gluc[keep]
 ins = ins[keep]
 10000/( sqrt(gluc[1]*ins[1]*sum(gluc*c(1,2,2,2,1))*
   sum(ins*c(1,2,2,2,1))/64 ) )
}
#matsudaPure ( gluc= c(91 ,130 ,148 ,157 ,157),
#   ins=c(6.9,24.0,37.6,59.1,71.3), times=c(0,30,60,90,120) )

#' Compute Matsuda's index (120 minute form) for all samples
#' @param mae instance of \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#' @param gname name of ExperimentList component holding glucose concentrations
#' @param iname name of ExperimentList component holding insulin concentrations
#' @param outname name of ExperimentList component holding Matsuda index
#' @param allowGaps logical, should addAssay fail if there are gaps between sample labels for new assay and rownames of colData
#' @return instance of \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}} that includes assay \code{"Mats120"}
#' @note The formula of the WEB CALCULATOR at \url{http://mmatsuda.diabetes-smc.jp/english.html} is used
#' @export
addMatsuda120 = function(mae, gname = "glucose", iname="insulin", outname="Mats120", allowGaps=FALSE) {
  stopifnot(class(mae) == "ogttCohort")
  times120 = c(0,30,60,90,120)
  stopifnot(all(times120 %in% mae@times))
  gluc = experiments(mae)[[gname]]
  ins = experiments(mae)[[iname]]
  times = mae@times
  stopifnot(length(times)>0)
  mats = matrix(sapply(1:ncol(gluc), function(x)
     matsuda120( gluc[,x], ins[,x], times )),nrow=1)
  rownames(mats) = outname
  colnames(mats) = colnames(gluc)
#  experiments(mae)[[outname]] = mats
#  mae
  addAssay(mae, mats, outname, allowGaps=allowGaps) # properly deals with sampleMap
}
