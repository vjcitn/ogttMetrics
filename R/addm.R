#' update an ogttCohort with insulin sensitivity by minimal model, with convergence information
#' @param mae instance of \code{\link{ogttCohort-class}}
#' @param iter lapply-like function, possibly mclapply
#' @param outSIname character token to name the assay in output MAE
#' @param replace logical, if FALSE will fail if "SI" present among names of assays(mae)
#' @param allowGaps logical as used in \code{\link{addAssay}}
#' @param \dots passed to \code{\link{getMinmodSIs}}
#' @examples
#' data(obaSamp)
#' lk = addMinmodSIs(obaSamp[,1:5])
#' lk
#' @export
 addMinmodSIs = function(mae, iter = lapply, outSIname="SI", replace=TRUE, 
   allowGaps=FALSE, ...) {
   nm = names(experiments(mae))
   if ("SI" %in% nm) {
     if (!replace) stop("replace is FALSE, but SI found among assays(mae)")
     mae = subsetByAssay(mae, setdiff(nm, "SI"))
     warning("SI found in mae, removing and replacing")
     }
   mat = getMinmodSIs(mae, iter, ...)
   addAssay(mae, mat, outSIname, allowGaps=allowGaps)
 }

.homa = function(glu, ins) {
 22.5 * 18 / ( glu[1] * ins[1] )
}
.quicki = function(glu, ins) {
 1 / ( log(glu[1]) + log(ins[1]) )
}

# any scalar index based on vectors glu and ins can be added to
# ogttCohort using this function
#
getScalar = function(oc, index=.homa, name="HOMA") {
 stopifnot(class(oc) == "ogttCohort")
 stopifnot(oc@times[1] == 0)
 g = assays(oc)$glucose
 ins = assays(oc)$insulin
 h = sapply(1:ncol(g), function(x) index(g[,x], ins[,x]))
 h = matrix(h, nrow=1)
 colnames(h) = colnames(oc)[[1]]
 rownames(h) = name
 h
}

#' update an ogttCohort with HOMA as defined in 10.1210/jc.2002-021127
#' @param oc ogttCohort instance
#' @examples
#' data(obaSamp)
#' addHOMA(obaSamp)
#' @export
addHOMA = function(oc) {
 addAssay(oc, getScalar(oc, .homa, "HOMA"), "HOMA")
}
 
#' update an ogttCohort with QUICKI as defined in 10.1210/jc.2002-021127
#' @param oc ogttCohort instance
#' @examples
#' data(obaSamp)
#' addQUICKI(obaSamp)
#' @export
addQUICKI = function(oc) {
 addAssay(oc, getScalar(oc, .quicki, "QUICKI"), "QUICKI")
}
