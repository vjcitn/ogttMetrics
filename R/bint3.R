## based on xoverSamp, add diet-specific Matsuda
#' @import magrittr
#' @importFrom dplyr inner_join
updateXMats120 = function(xoverMAE) {
  stopifnot(length(xoverMAE@times)>0)
  xoverMAE = addMatsuda120(xoverMAE, 
      "glucose_CG", "insulin_CG", "Mats120_CG", allowGaps=TRUE)
  xoverMAE = addMatsuda120(xoverMAE, 
      "glucose_Cg", "insulin_Cg", "Mats120_Cg", allowGaps=TRUE)
  xoverMAE = addMatsuda120(xoverMAE, "glucose_cG", 
      "insulin_cG", "Mats120_cG", allowGaps=TRUE)
  xoverMAE = addMatsuda120(xoverMAE, 
      "glucose_cg", "insulin_cg", "Mats120_cg", allowGaps=TRUE)
  xoverMAE
}

updateXSI = function(xoverMAE) {
  stopifnot(length(xoverMAE@times)>0)
  xoverMAE = addMinmodSIs(xoverMAE, 
      gname="glucose_CG", iname="insulin_CG", 
      outSIname="SI_CG", allowGaps=TRUE)
  xoverMAE = addMinmodSIs(xoverMAE, 
      gname="glucose_Cg", iname="insulin_Cg", 
      outSIname="SI_Cg", allowGaps=TRUE)
  xoverMAE = addMinmodSIs(xoverMAE, 
      gname="glucose_cG", iname="insulin_cG", 
      outSIname="SI_cG", allowGaps=TRUE)
  xoverMAE = addMinmodSIs(xoverMAE, 
      gname="glucose_cg", iname="insulin_cg", 
      outSIname="SI_cg", allowGaps=TRUE)
  xoverMAE
}

## code diet and rowname in xoverSamp, after updateXMats120
harmonizeX = function(xoverMAE, stub="Mats120") {
  ll = longFormat(xoverMAE)
  ll$rowname = gsub(paste0(stub, "_.*"), stub, ll$rowname)
  ll$diet = gsub( ".*_(.*)"  , "\\1"  ,  as.character(ll$assay))
  ll
}
  
## merge two diets in an updated xoverMAE in long form
  merge_diets = function(lth, d1="CG", d2="cg", measure="Mats120") {
    di1 = data.frame(lth) %>% 
           dplyr::filter(diet==d1 & rowname==measure) %>% 
           dplyr::rename(ocd1 = value)
    di2 = data.frame(lth) %>% 
           dplyr::filter(diet==d2 & rowname==measure) %>% 
           dplyr::rename(ocd2 = value)
    inner_join(di1, di2, by="primary")
  }
  
## given a long form updated xoverMAE, generate pairwise diet contrasts
generateTests = function(lth, pairs, levels, measure="Mats120") {
    allm = lapply(pairs, function(x) merge_diets(lth, d1=x[1], d2=x[2], measure=measure))
    ans = lapply(seq_along(allm), function(x) with(allm[[x]], t.test(ocd1-ocd2,
     var.equal=TRUE, conf.level=levels[x])))
    nms = sapply(pairs, paste, collapse="-")
    names(ans) = nms 
    ans 
  }

#' generate statistics for analog of OMNICarb paper figure 3
#' @param xoverMAE MAE instance like xoverSamp in data
#' @param contrasts a list of diet code pairs
#' @param levels confidence coefficients to be used
#' @param type character either "Mats120" or "SI"
#' @examples
#' data(xoverSamp)
#' tt = fig3tests(xoverSamp)
#' fig3plot(tt)
#' @export
fig3tests = function(xoverMAE, 
   contrasts= list(c("Cg", "CG"), c("cg", "cG"), c("cG", "CG"), 
      c("cg", "Cg"), c("cg", "CG")), 
   levels=c(.95, .95, .95, .95, .99), type="Mats120") {
  if (type=="Mats120")
       xoverMAE = updateXMats120(xoverMAE)
  else if (type=="SI")
       xoverMAE = updateXSI(xoverMAE)
  ll = harmonizeX(xoverMAE, stub=type)
  curtests = generateTests(ll, contrasts, levels, measure=type )
  ints = sapply(curtests, "[[", "conf.int")
  mns = sapply(curtests, "[[", "estimate")
  nint = length(mns)
  list(nint=nint, tests=curtests, ints=ints, mns=mns)
}

#' plot analog of OMNICarb paper figure 3
#' @param f3out output of fig3tests
#' @param measTag character tag to be used in xlab after delta symbol
#' @export
fig3plot = function (f3out, measTag = "Matsuda") {
    xl = range(as.numeric(f3out$ints)) * 1.05
    par(mar = c(4, 5, 2, 2))
    if (measTag == "Matsuda") xlb = expression(paste(Delta, " Matsuda"))
    else if (measTag == "SI") xlb = expression(paste(Delta, " SI"))
    plot(f3out$mns, f3out$nint:1, pch = 19, axes = FALSE, xlab=xlb,
        xlim = xl, ylim = c(0.8, f3out$nint + 
        0.2), ylab = " ")
    segments(f3out$ints[1, ], f3out$nint:1, f3out$ints[2, ], 
        f3out$nint:1)
    axis(1)
    axis(2, at = f3out$nint:1, labels = names(f3out$tests), las = 2)
    abline(v = 0)
}

fig3plotOLD = function(f3out, measTag = "Matsuda") {
  xl = range(as.numeric(f3out$ints))*1.05
  par(mar=c(4,5,2,2))
  plot(f3out$mns, f3out$nint:1, pch=19, axes=FALSE, 
      xlab=expression(paste(Delta, paste0(" ",measTag))), xlim=xl,
      ylim=c(.8,f3out$nint+.2), ylab=" ")
  segments(f3out$ints[1,], f3out$nint:1, f3out$ints[2,], f3out$nint:1)
  axis(1)
  axis(2, at=f3out$nint:1, labels=names(f3out$tests),las=2)
  abline(v=0)
}
  
