#' combine selected elements of a group of MultiAssayExperiments
#' @param maes a list of MultiAssayExperiment instances
#' @param assays2keep a character vector, used with subsetByAssay to restrict the elements of \code{maes}
#' @param basenames a character vector giving generic names for assay elements; if null, uses names of ElementList(maes)[[1]]
#' @param \dots not currently used
#' @export
concat = function( maes = list(
  base=omnicBase, 
  CG=omnicCG, Cg=omnicCg, Cg=omnicCg, cg=omniccg
  ), assays2keep = NULL, basenames=NULL, ...) {
 pd = pData(maes[[1]])
 if (is.null(assays2keep))
     exs = lapply(maes, function(x) experiments(x))
    else exs = lapply(maes, function(x) experiments(subsetByAssay(x, assays2keep)))
 numel = sum(sapply(exs, length))
 enames = names(maes)
 ans = vector("list", numel) # flat list
 k = 1
 for (i in 1:length(exs)) {
   for (j in 1:length(exs[[i]])) {
     ans[[k]] = exs[[i]][[j]]
#     names(ans[[k]]) = paste(names(exs[[i]][[j]]), enames[i], sep="_")
     k = k+1
     }
   }
 if (is.null(basenames)) basenames = names(exs[[1]])
 assn = expand.grid(basenames, enames)
 names(ans) = apply(t(assn),2,paste, collapse="_")
# FIXME -- should inherit
 MultiAssayExperiment( ExperimentList(ans), pData=pd )
}

generateTests = function(lth, pairs, levels, measure="Mats120",
  testDiffs=function(delta,...) t.test(delta, var.equal=TRUE, conf.level=.95), ...) {
  allm = lapply(pairs, function(x) 
     merge_diets(lth, d1=x[1], d2=x[2], measure=measure))
  ans = lapply(seq_along(allm), function(x) with(allm[[x]], testDiffs(ocd1-ocd2,
   ...)))
  nms = sapply(pairs, paste, collapse="-")
  names(ans) = nms
  ans
}

