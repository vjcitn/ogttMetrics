#' @import methods
#' @import knitr
#' @import MultiAssayExperiment
#' @import S4Vectors
#' @import Biobase
#' @import SummarizedExperiment
#' @import deSolve
#' @importFrom stats approxfun
#' @importFrom stats coef
#' @importFrom stats nls
#' @importFrom stats nls.control
#' @importFrom stats prcomp
#' @importFrom stats predict
#' @importFrom stats t.test
#' @importFrom utils sessionInfo
#' @importFrom graphics par
#' @importFrom graphics lines
#' @importFrom graphics plot
#'

#'
#' @rdname ogttCohort
#' @title 
#' an extension of MultiAssayExperiment including timing information
#' @name ogttCohort-class
#' @slot times a numeric vector of OGTT sampling times
#' @exportClass ogttCohort
setClass("ogttCohort", representation(times="numeric"),
   contains="MultiAssayExperiment")
#
setAs("MultiAssayExperiment", "ogttCohort", function(from) {
  stopifnot(!is.null(metadata(from)$times))
  new("ogttCohort", times=metadata(from)$times,
    ExperimentList=from@ExperimentList,
    pData=from@pData,
    sampleMap=from@sampleMap,
    metadata=from@metadata,
    drops=from@drops)
})
#'
#' add a matrix-valued assay to a MultiAssayExperiment instance
#' @param mae MultiAssayExperiment instance
#' @param mat a matrix with rownames assigned; if colnames absent will be taken from rownames(pData(mae))
#' @param assayname character string used to name assay in ExperimentList of mae
#' @examples
#' data(obaSamp)
#' lit = subsetByAssay(obaSamp, c("glucose", "insulin"))
#' d = dim(ee <- experiments(lit)$glucose)
#' dum = matrix(rnorm(prod(d)), nrow=d[1])
#' rownames(dum) = rownames(ee)
#' colnames(dum) = colnames(ee)
#' addAssay(lit, dum, "dum")
#' @export
addAssay = function(mae, mat, assayname) {
 stopifnot(class(mae)=="ogttCohort") # not yet available for MAE
 stopifnot(!missing(assayname))
 stopifnot(is(rownames(mat), "character"))
 if (is.null(colnames(mat))) {
   warning("assigning colnames to mat")
   colnames(mat) = rownames(pData(mae))
 } else {
   stopifnot(all.equal(colnames(mat), rownames(pData(mae))))
   }
 nl = list(mat)
 names(nl) = assayname
 newEL = ExperimentList(c(experiments(mae)@listData, nl))
 inmd = metadata(mae)  # will be NULL if mae is ogttCohort
 newmae = MultiAssayExperiment(newEL, pData(mae)) # set up sampleMap, lose ogttCohort aspect
 new("ogttCohort", newmae, times=mae@times)
}
# metadata(newmae) = inmd
# if (is(mae, "ogttCohort")) {  # restore ogttCohort aspect (times slot)
#   newmae = as(newmae, "ogttCohort")
#   newmae@times = mae@times
#   }
# newmae
#}