# for R 3.3
rearrangeEL = function (object, shape = "long", ...) 
{
    dataList <- as.list(object)
    dataList <- lapply(seq_along(object), function(i, flatBox) {
        S4Vectors::DataFrame(assay = S4Vectors::Rle(names(object)[i]), 
            rearrange(flatBox[[i]], ...))
    }, flatBox = object)
    dataList
}

rearrangeMAE = function (object, shape = "long", pDataCols = NULL, 
        ...) 
    {
        addCols <- !is.null(pDataCols)
        dataList <- rearrangeEL(experiments(object), ...)
        dataList <- lapply(dataList, function(rectangleDF) {
            primary <- S4Vectors::Rle(sampleMap(object)[match(rectangleDF[["colname"]], 
                sampleMap(object)[["colname"]]), "primary"])
            rectangleDF <- S4Vectors::DataFrame(rectangleDF, 
                primary = primary)
            rectangleDF[, c("assay", "primary", "rowname", "colname", 
                "value")]
        })
        outputDataFrame <- do.call(rbind, dataList)
        if (addCols) {
            extraColumns <- pData(object)[, pDataCols, drop = FALSE]
            rowNameValues <- rownames(extraColumns)
            rownames(extraColumns) <- NULL
            matchIdx <- BiocGenerics::match(outputDataFrame[["primary"]], 
                rowNameValues)
            outputDataFrame <- BiocGenerics::cbind(outputDataFrame, 
                extraColumns[matchIdx, , drop = FALSE])
        }
        if (shape == "wide") {
            outputDataFrame <- as.data.frame(outputDataFrame)
            outputDataFrame <- tidyr::unite_(outputDataFrame, 
                "feature", c("assay", "rowname", "colname"))
            outputDataFrame <- tidyr::spread(outputDataFrame, 
                key = "feature", value = "value")
            outputDataFrame <- DataFrame(outputDataFrame)
        }
        return(outputDataFrame)
    }

#' get multivariate outlier indices separately for glucose and insulin series
#' @importFrom parody mv.calout.detect
#' @param oc ogttCohort instance
#' @param t_glu function to transform glucose data to approximate multivariate normality
#' @param t_ins function to transform insulin data to approximate multivariate normality
#' @param \dots passed to \code{\link[parody]{mv.calout.detect}}
#' @examples
#' data(obaSamp)
#' mvOutliers(obaSamp)
#' mvOutliers(obaSamp, alpha=.1)
#' @export
mvOutliers = function(oc, t_glu=force, t_ins=force, ...) {
 requireNamespace("parody")
 a = assay(oc)
 ins = na.omit(t(a$insulin))
 if (!is.null(nd <- attributes(ins)$na.action)) warning(paste0(length(nd), " records dropped with NA in insulin"))
 glu = na.omit(t(a$glucose))
 if (!is.null(nd <- attributes(glu)$na.action)) warning(paste0(nd, " records dropped with NA in glucose"))
 insout = mv.calout.detect(t_ins(ins))
 gluout = mv.calout.detect(t_glu(glu))
 ids = colnames(oc)
 insids = ids$insulin
 gluids = ids$glucose
 if (!is.na(insout$ind[1])) insout$ind = insids[insout$ind]
 if (!is.na(gluout$ind[1])) gluout$ind = gluids[gluout$ind]
 list(insulinOutliers=insout$ind, glucoseOutliers=gluout$ind)
}

ogbox = function(oc, type="glucose") {
 thin = rearrangeMAE(oc)
 thin$assay = as.character(thin$assay)
 thin = as.data.frame(thin[which(thin$assay==type),])
 if (type=="glucose") thin$rowname = gsub("gluc", "", gsub("gluct", "", thin$rowname))
 else if (type=="insulin") thin$rowname = gsub("ins", "", thin$rowname)
 levs = oc@times
 thin$rowname = factor(thin$rowname, levels=as.character(levs))
 ggplot(thin, aes(x=rowname, y=value)) + geom_boxplot() + ylab(type) + 
   xlab("min") + theme_gray()
}

#' simple quality assessments for OGTT series in a cohort
#' @importFrom ggbiplot ggbiplot
#' @param oc ogttCohort instance
#' @param choices a 2-vector of integers passed as \code{choices} to \code{\link[ggbiplot]{ggbiplot}}, selecting principal axes for biplot display
#' @examples
#' example(csvImport) # makes democ
#' QCplots(democ)
#' @export
QCplots = function(oc, choices=1:2) {
 requireNamespace("parody")
 requireNamespace("cowplot")
 requireNamespace("ggplot2")
 a = assay(oc)
 ins = na.omit(data.frame(t(a$insulin)))
 if (!is.null(nd <- attributes(ins)$na.action)) warning(paste0(length(nd), " records dropped with NA in insulin"))
 inslab = colnames(oc)$insulin
 if (length(nd)>0) inslab=inslab[-nd]
 glu = na.omit(data.frame(t(a$glucose)))
 if (!is.null(nd <- attributes(glu)$na.action)) warning(paste0(length(nd), " records dropped with NA in glucose"))
 glulab = colnames(oc)$glucose
 if (length(nd)>0) glulab=glulab[-nd]
 pins = prcomp(ins)
 pglu = prcomp(glu)
 CH1 = gsub("%%N%%", choices[1], "PC%%N%%")
 CH2 = gsub("%%N%%", choices[2], "PC%%N%%")
 glubi = ggbiplot(pglu,choices=choices,labels=glulab) + xlab(CH1) + ylab(CH2) + theme_gray()
 insbi = ggbiplot(pins,choices=choices,labels=inslab) + xlab(CH1) + ylab(CH2)+ theme_gray()
 plot_grid(ogbox(oc), ogbox(oc, "insulin"), glubi, insbi, nrow=2, ncol=2) 
}
