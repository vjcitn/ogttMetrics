#' support for importing separate CSV with fat format: one record per individual
#' @importFrom S4Vectors DataFrame
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList
#' @param gpath character pathname of glucose data
#' @param ipath character pathname of insulin data
#' @param spath character pathname of sample-level data
#' @param times numeric vector of observation times in minutes
#' @param gfilt function that operates on glucose data accepting and returning data.frame as returned by read.csv
#' @param ifilt function that operates on insulin data accepting and returning data.frame as returned by read.csv
#' @param sfilt function that operates on sample data accepting and returning data.frame as returned by read.csv
#' @examples
#' pref = system.file("csv_example", package="ogttMetrics")
#' gpath = paste0(pref, "/glucBase.csv")
#' ipath = paste0(pref, "/insBase.csv")
#' spath = paste0(pref, "/sampBase.csv")
#' democ = csvImport(gpath, ipath, spath)
#' democ
#' @export
csvImport = function(gpath, ipath, spath, times=c(0,10,20,30,60,90,120),
  gfilt=force, ifilt=force, sfilt=force) {
 #
#"id","gluc0","gluc10","gluc20","gluc30","gluc60","gluc90","gluc120"
#
 gluc = gfilt(read.csv(gpath))
 ins = ifilt(read.csv(ipath))
 samps = DataFrame(sfilt(read.csv(spath)))
 glucCol1Name = names(gluc)[1]
 glucCol2Name = names(gluc)[2]
 insCol1Name = names(ins)[1]
 sampCol1Name = names(samps)[1]
 stopifnot(glucCol1Name=="id", 
        insCol1Name=="id",
        sampCol1Name=="id",
        glucCol2Name=="gluc0",
        ncol(gluc)==length(times)+1, 
        ncol(ins)==length(times)+1,
        nrow(gluc)==nrow(ins),
        nrow(gluc)==nrow(samps),
        all.equal(gluc$id, ins$id),
        all.equal(gluc$id, samps$id))
 colntimes = as.integer(gsub("gluc", "", names(gluc)[-1]))
 stopifnot(all.equal(colntimes, as.integer(times)))
 glucose=data.matrix(t(gluc[,-1]))
 insulin=data.matrix(t(ins[,-1]))
 colnames(glucose) = gluc$id
 colnames(insulin) = gluc$id
 rownames(samps) = samps$id
 el = ExperimentList(list(glucose=glucose, insulin=insulin))
 tmp = MultiAssayExperiment(el, colData=samps)
 new("ogttCohort", tmp, times=times)
}

csvImport2 = function(gpath, ipath, spath, times=c(0,10,20,30,60,90,120),
  gfilt=force, ifilt=force, sfilt=force) {
 #
#"id","gluc0","gluc10","gluc20","gluc30","gluc60","gluc90","gluc120"
#
 gluc = gfilt(read.csv(gpath))
 ins = ifilt(read.csv(ipath))
 samps = DataFrame(sfilt(read.csv(spath)))
 glucCol1Name = names(gluc)[1]
 insCol1Name = names(ins)[1]
 sampCol1Name = names(samps)[1]
 stopifnot(glucCol1Name=="id", 
        insCol1Name=="id",
        sampCol1Name=="id",
        ncol(gluc)==length(times)+1, 
        ncol(ins)==length(times)+1,
        nrow(gluc)==nrow(ins),
        nrow(gluc)==nrow(samps),
        all.equal(gluc$id, ins$id),
        all.equal(gluc$id, samps$id))
 glucose=data.matrix(t(gluc[,-1]))
 insulin=data.matrix(t(ins[,-1]))
 list(gluc=glucose, ins=insulin, samps=samps)
}
