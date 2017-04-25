geta6=function(a1,a2,a3,BW,D=75000,FA=.87,DC=.017){
        area1=(a1*30)/2
        area2=((a1+a2)*30)/2
        area3=((a2+a3)*30)/2
        area90=area1+area2+area3
        area=(D*FA)/BW
        missing=area-area90
        e420=(exp(-DC*300)/-DC)
        e120=(exp(-DC*0)/-DC)
        a6=(missing-15*a3)/(15+e420-e120)
        a6;
        }

ra=function(a1,a2,a3,t,BW,D=75000,FA=.87,DC=.017){
        a6=geta6(a1,a2,a3,BW,D,FA,DC)

        if(t<30){ra=((a1/30)*t)}
        if((t>=30)&(t<60)){ra=(a1+((a2-a1)/30)*(t-30))}
        if((t>=60)&(t<90)){ra=(a2+((a3-a2)/30)*(t-60)) }
        if((t>=90)&(t<120)){ra=(a3+((a6-a3)/30)*(t-90)) }
        if(t>=120){ra=(a6*exp(-DC*(t-120)))}
        ra
}


#' fit the dalla Man version of the Bergman minimal model (piecewise linear form for rate of appearance of glucose)
#' @param insulin a numeric vector of insulin values recorded at times in \code{times}
#' @param glucose a numeric vector of glucose concentrations, same length as \code{insulin}
#' @param bodyweight a numeric scalar
#' @param times times in minutes
#' @param nlsmaxit maximum number of iterations for NLS
#' @param nlstol tolerance for convergence of NLS
#' @param fullNLScontrol an optional list of control parameters for NLS
#' @param nlstrace logical requesting retention of trace
#' @param Sg fractional glucose effectiveness (default 0.028 min^-1)
#' @param V volume of distribution (default 1.34 dl/kg)
#' @param D ingested glucose dose in mg (default 75000)
#' @param FA fraction of glucose dose that is absorbed, default = 0.87, dimensionless
#' @param DC decay constant for model for rate of appearance of glucuse, default 0.017
#' @param p2 rate constant of the remote insulin compartment from which insulin action emanates, default value 0.012 min^-1
#' @param nlsinit initial values for unknown parameters a1, a2, a3 (coefficients of piecewise linear model for rate of appearance of glucose) and SI (insulin sensitivity)
#' @param trajTimes a vector of times at which glucose and insulin trajectories should be estimated
#' @param inputID a string that may be used for rendering an ID tag
#' @return a list with elements \code{fit} for the NLS fit, \code{soln} for the lsoda based solution to the differential equation model, \code{pwlinRa} for the piecewise-linear estimate of Ra (rate of appearance of glucose), \code{input} for a sublist with information about the input data, and \code{sessInf} for the sessionInfo metadata
#' @export
fitOneMinMod = function (insulin, glucose, bodyweight, times, 
   nlsmaxit = 100, nlstol=.1, 
   fullNLScontrol=NULL, nlstrace=FALSE,
   Sg=0.028, V=1.34, D=75000, FA=.87, DC=.017, p2=0.012, 
   nlsinit=list(a1=4, a2=4, a3=4, SI=1e-4), 
   trajTimes=0:120, inputID="") {
#
# code originated by BJ Harshfield, some enhancements by VJC
#
    basicNLScontrol = list(warnOnly=TRUE,
         maxiter = nlsmaxit, tol = nlstol, minFactor=1e-12)
    if (is.null(fullNLScontrol)) fullNLScontrol = basicNLScontrol
     else stopifnot(all(c("warnOnly", "maxiter", "tol", "minFactor") %in%
                names(fullNLScontrol)))
    BW = bodyweight
    g = glucose
    nag = which(is.na(g))
    i = insulin
    nai = which(is.na(i))
    t = times
    if (length(c(nag, nai))>0) {
          drop = sort(unique(c(nag, nai)))
          g = g[-drop]
          i = i[-drop]
          t = t[-drop]
          message("missing values found in assay series, omitted")
          }
    Gb = as.numeric(g[1])
    Ib = as.numeric(i[1])
    Insulin = approxfun(t, i, rule = 2)
    model <- function(t, Y, parameters) {
        with(as.list(parameters), {
            dy1 = -(Sg + Y[2]) * Y[1] + Sg * Gb + (ra(a1, a2, 
                a3, t, BW, D, FA, DC)/V)
            dy2 = -p2 * Y[2] + p2 * SI * (Insulin(t) - Ib)
            list(c(dy1, dy2))
        })
    }
    mmsolfn = function(a1, a2, a3, SI) 
        lsoda(c(Gb, 0), t, model, c(a1 = a1, a2 = a2, 
          a3 = a3, SI = SI))
    fit <- nls(g ~ mmsolfn(a1, a2, a3, SI)[,2], 
        start = nlsinit,
        trace = nlstrace, control = fullNLScontrol)
    co = coef(fit)
    soln = lsoda(c(Gb, 0), trajTimes, model, co)
    colnames(soln) = c("t", "G(t)", "X(t)")
    pwlinRa = sapply(trajTimes,
      function(curt) ra(co[1],co[2],co[3],curt,BW=BW,D=D,FA=FA,
         DC=DC))
    ans = list(fit=fit, soln=soln, pwlinRa = pwlinRa,
       input=list(glucose=glucose, insulin=insulin, times=times, bodyweight=bodyweight, trajTimes=trajTimes), sessInf=sessionInfo(), ID=inputID)
    class(ans) = "oneMinMod"
    ans
}

#' @rdname minmodByID
#' @param x instance of S3 class oneMinMod
#' @export
print.oneMinMod = function(x, ...) {
  cat("ogttMetrics Minimal Model fit:\n")
  print(x$fit)
  cat("---\n  other components are: ", names(x)[-1], "\n")
}
#'
#' Fit minimal model for a single subject
#' @param mae MultiAssayExperiment instance
#' @param id character id found in rownames(colData(mae))
#' @param gname name of ExperimentList component holding glucose concentrations
#' @param iname name of ExperimentList component holding insulin concentrations
#' @param \dots passed to \code{\link{fitOneMinMod}}
#' @examples
#' data(obaSamp)
#' m1 = minmodByID(obaSamp, "1")
#' m1
#' @export
minmodByID = function(mae, id, gname="glucose", iname="insulin",  ...) {
 id = as.character(id)
 stopifnot(id %in% rownames(colData(mae)))
 times = mae@times #metadata(mae)$times
 ae = MultiAssayExperiment::experiments
 gluc = ae(subsetByAssay(mae, gname))[[1]][,id]
 stopifnot(length(times) == length(gluc))
 ins = ae(subsetByAssay(mae, iname))[[1]][,id]
 bw = ae(subsetByAssay(mae, "bodyweight"))[[1]][,id]
 ans = fitOneMinMod( insulin=ins, glucose=gluc, 
          bodyweight=bw, times=times, ... )
 ans$input$ID = id
 ans$ID = id # need to prune
 class(ans) = "minmodByID"
 ans
}

#' @rdname minmodByID
#' @export
print.minmodByID = function(x, ...) {
 cat("minmodByID result for ID ", x$input$ID, "\n", sep="")
 cat("SI estimate: ", coef(x$fit)["SI"], "\n")
 cat("object has components ", names(x))
 cat("\n")
}

getMinmodSIs_obsolete = function(mae, iter=lapply) {
 allid = rownames(colData(mae))
 ans = iter(allid, function(x) { try(minmodByID( mae, x))})
 sapply(ans, function(x)
   {
   if(inherits(x, "try-error")) return(NA) else coef(x$fit)["SI"]
   })
}

#'
#' Obtain SI from fits for all participants in an OGTT MultiAssayExperiment
#' @param mae MultiAssayExperiment instance
#' @param iter an lapply analog, supply mclapply for multicore (with options(mc.cores=...) appropriately set)
#' @param \dots passed to \code{\link{minmodByID}}
#' @export
getMinmodSIs = function (mae, iter = lapply, ...) 
{
    allid = rownames(colData(mae))
    ans = iter(allid, function(x) {
        try(minmodByID(mae, x, ...))
    })
    ans = sapply(ans, function(x) {
        if (inherits(x, "try-error")) 
            return(rep(NA,5))
        else {
             convdat = x$fit$convInfo[1:4] # avoid 'list' character owing to mix of logical and numeric
             nc = names(convdat)
             convdat = as.numeric(convdat)
             names(convdat) = nc
	     c(SI=as.numeric(coef(x$fit)["SI"]), convdat)
             } 
    })
    colnames(ans) = colnames(mae)$glucose
    ans
}

fmm4pop = function (insulin, glucose, bodyweight, times, id, SI,
   nlsmaxit = 100, nlstol=.1, 
   fullNLScontrol=nls.control(maxiter=3, warnOnly=TRUE), nlstrace=TRUE,
   Sg=0.028, V=1.34, D=75000, FA=.87, DC=.017, p2=0.012, 
   nlsinit=list(a1=4, a2=4, a3=4, SI=SI), 
   trajTimes=0:120, inputID="") {
#
# code originated by BJ Harshfield, some enhancements by VJC
#
    inslist = split(insulin, id)
    gluclist = split(glucose, id)
    bwlist = split(bodyweight, id)
    timeslist = split(times, id)
    nsubject = length(unique(id))
    SI = SI[1]
    unlist(lapply(1:nsubject, function(ind) {
       predict(fitOneMinMod(inslist[[ind]], glucose[[ind]], bwlist[[ind]][1],
             timeslist[[ind]], nlsinit=list(a1=4, a2=4, a3=4, SI=SI),
             fullNLScontrol=fullNLScontrol)$fit)
       }))
}
     
