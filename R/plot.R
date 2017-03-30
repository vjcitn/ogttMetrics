#' crude visualization developed prior to plot_OGTT_fit
#' @param x ogtt obj
#' @param \dots not used
#' @export
#' @export
plot2 = function(x, ...) {
 opar = par(no.readonly=TRUE)
 on.exit(par(opar))
 par(mfrow=c(2,2), mar=c(4,4,2,1))
 tim = x$input$times
 glu = x$input$glucose
 ins = x$input$insulin
 nag = which(is.na(glu))
 nai = which(is.na(ins))
 if (length(c(nag,nai))>0) {
    drop = sort(unique(c(nag,nai)))
    tim = tim[-drop]
    ins = ins[-drop]
    glu = glu[-drop]
    }
 plot(tim, glu, main=paste0("ID ", x$input$ID),
   ylim=c(min(glu-5), max(x$soln[,"G(t)"])+5),
   xlab="t (min)", ylab="G(t)")
 lines(x$input$trajTimes, x$soln[, "G(t)"])
 plot(tim, ins, xlab="t (min)", ylab="I(t)")
 plot(x$soln[,"X(t)"], x$soln[, "G(t)"], type="l", lty=1,
   ylim=c(min(x$soln[,"G(t)"])-5, max(x$soln[,"G(t)"])+5),
   xlim=c(min(x$soln[,"X(t)"])-.01, max(x$soln[,"X(t)"])+.01),
   xlab="X(t)", ylab="G(t)")
 plot(x$input$trajTimes, x$pwlinRa, type="l", main="Ra(t)", xlab="t", ylab="Ra(t)")
 invisible(NULL)
}

#' produce a 4-panel plot of aspects of the OGTT minimal model fit
#' @param x output of minmodByID
#' @param ptitles a list of strings, with names 'a', ..., 'd', that will
#'  be used as title elements for each panel
#' @param \dots not used
#' @return northwest panel is a scatter plot of observed glucose vs. time,
#' with a solid line sketching the trajectory of the fitted mean glucose
#' based on the NLS fit; northeast panel is a scatterplot of observed insulin
#' vs. time; southwest panel is a sketch of fitted mean glucose against the
#' abstract variable X(t) [insulin action on glucose production and disposal];
#' southeast panel is a sketch of the fitted piecewise linear model for Ra(t),
#' the rate of appearance of glucose in the OGTT
#' @importFrom ggplot2 ggplot geom_point aes geom_line ggtitle xlab ylab theme_gray geom_text geom_label
#' @importFrom cowplot plot_grid
#' @examples
#' data(obaSamp)
#' m1 = minmodByID(obaSamp, "1")
#' plot_OGTT_fit(m1)
#' @export
plot_OGTT_fit = function(x, 
   ptitles=list(a="(a)", b="(b)", c="(c)", d="(d)"), ...) {
 tim = x$input$times
 glu = x$input$glucose
 ins = x$input$insulin
 nag = which(is.na(glu))
 nai = which(is.na(ins))
 if (length(c(nag,nai))>0) {
    drop = sort(unique(c(nag,nai)))
    tim = tim[-drop]
    ins = ins[-drop]
    glu = glu[-drop]
    }
 df = data.frame(time=x$input$trajTimes, Gt=x$soln[,"G(t)"])
 pt_df = data.frame(ins=ins, glu=glu, tim=tim)
 g1 = ggplot(df, aes(x=time, y=Gt)) + geom_line() + 
        geom_point(aes(x=tim, y=glu), data=pt_df) + xlab("t (min)") +
        ylab("G(t)") + theme_gray() + ggtitle(ptitles[["a"]])
 g2 = ggplot(pt_df, aes(x=tim, y=ins)) + geom_line() + geom_point() +
        xlab("t (min)") + ylab("I(t)") + theme_gray() + ggtitle(ptitles[["b"]])
 dy_df = data.frame(Xt=x$soln[,"X(t)"], Gt=x$soln[, "G(t)"])
 ld = dy_df[as.integer(tim)+1,]
 rownames(ld) = as.character(as.integer(tim))
 ld$col = factor(as.character(tim), levels=as.character(tim))
 ld$min = factor(as.character(tim), levels=as.character(tim))
 ld = ld[c("0", "60", "120"),]
 g3 = ggplot(dy_df, aes(y=Gt, x=Xt)) + geom_point() +
        ylab("G(t)") + xlab("X(t)") + 
        theme_gray() + ggtitle(ptitles[["c"]]) +
        geom_label(data = ld, aes(x = Xt, y = Gt, 
        label = rownames(ld), colour=min, fill=min), size = 1)
 radf = data.frame(ra=x$pwlinRa, t=x$input$trajTimes)
 g4 = ggplot(radf, aes(y=ra, x=t)) + geom_line() + xlab("t (min)") +
    ylab("Ra(t)") + theme_gray() + ggtitle(ptitles[["d"]])
 plot_grid(g1, g2, g3, g4, nrow=2)
}
