
 dumpOne = function(oc, id, ...) {
  g = experiments(oc)$glucose[,id]
  ins = experiments(oc)$insulin[,id]
  bw = experiments(oc)$bodyweight[,id]
  writeLines("DATA", ...)
  writeLines(paste("CONST weight ", round(as.numeric(bw),2), sep=""), ...)
  writeLines("T glucose insulin", ...)
  for (i in 1:length(g))
    writeLines(paste(c(oc@times[i], round(g[i], 3), round(ins[i], 3)), collapse=" "), ...)
  }

#' an ad hoc tool for preparing individual files for analysis in SAAM-II
#' @param oc instance of ogttCohort
#' @examples
#' data(obaSamp)
#' owd = getwd()
#' td = tempdir()
#' setwd(td)
#' dumpOc(obaSamp[,1:3])
#' readLines("ogtt_1.dat")
#' setwd(owd)
#' @export
dumpOc = function(oc) {
 ids = colnames(oc)[[1]]
 for (i in ids) {
   f = file(paste("ogtt_", i, ".dat", sep=""), "w")
   dumpOne( oc, i, f)
   close(f)
   }
}
