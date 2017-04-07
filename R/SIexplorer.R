#' interactive multipanel visualization for ogttCohort instance
#' @import shiny
#' @import ggbiplot
#' @importFrom plotly renderPlotly ggplotly plotlyOutput
#' @importFrom ggplot2 scale_x_sqrt scale_y_log10
#' @param oc ogttCohort instance
#' @param winsorizeSI if TRUE, move negative estimates of SI to smallest positive value
#' @param \dots passed to \code{\link{minmodByID}}
#' @examples
#' if (interactive()) {
#'   if (options()$example.ask) stop("must set options(example.ask=FALSE) before running example")
#'   data(obaSamp)
#'   SIexplorer(obaSamp)
#' }
#' @export
SIexplorer = function(oc=obaSamp, winsorizeSI=TRUE, ...) {
 stopifnot(is(oc, "ogttCohort"))
 eln = names(experiments(oc))
 stopifnot("Mats120" %in% eln)
 stopifnot("SI" %in% eln)
# get dataset name, sample names
 dstxt = deparse(substitute(oc))
 allids = colnames(oc)$glucose
#
# build up information for biplots
#
 a = assay(experiments(oc))
 ins = na.omit(data.frame(t(a$insulin)))
 insdrop = attributes(ins)$na.action
 insids = colnames(oc)[[1]]
 if (length(insdrop)>0) insids = colnames(oc)[[1]][-insdrop]
 glu = na.omit(data.frame(t(a$glucose)))
 gludrop = attributes(glu)$na.action
 gluids = colnames(oc)[[1]]
 if (length(gludrop)>0) gluids = colnames(oc)[[1]][-gludrop]
 pins = prcomp(ins)
 pglu = prcomp(glu)
#
# start app design
#
 ui = fluidPage(
       sidebarLayout(
        sidebarPanel(
         textOutput("datasetName"),
         textOutput("nrec"),
         helpText("Select ID for OGTT component plotting"),
         selectInput("idpick", "id", choices=allids,
                selected=allids[1]),
         selectInput("pcbipick_1", "PC for biplot X", choices=as.character(1:5),
                selected="1"),
         selectInput("pcbipick_2", "PC for biplot Y", choices=as.character(1:5),
                selected="2"),
         width=3), # end sidebarPanel
        mainPanel(
         tabsetPanel(
          tabPanel("SI vs Matsuda", 
            fluidRow(helpText("Cohort-wide SI vs Matsuda, hover over for ID, stats")),
            fluidRow(plotlyOutput("sivmat"))
             ),
          tabPanel("Glucose biplot", plotlyOutput("biplots_glucose")),
          tabPanel("Insulin biplot", plotlyOutput("biplots_insulin")),
          tabPanel("minmod fit", 
            fluidRow(textOutput("idtext")),
            fluidRow(textOutput("sitext")),
            fluidRow(plotOutput("demo"))
            ) # end fit panel
          )
         ) #end mainpanel
        ) # end layout
       ) # end page
 server = function(input, output, session) {
    output$datasetName = renderText(paste("Dataset:", dstxt))
    output$nrec = renderText(paste("# OGTT:", length(colnames(oc)[[1]])))
    output$idtext = renderText(paste("ID =", input$idpick))
    output$sitext = renderText(paste("SI =", round(assay(experiments(oc))$SI[1, input$idpick], 6)))
    output$sivmat = renderPlotly( {
       newdf = data.frame(id=colnames(oc)[[1]], 
                   mats120=as.numeric(assay(experiments(oc))$Mats120),
                   SI=as.numeric(assay(experiments(oc))$SI[1,]),
                   converged=as.logical(as.numeric(assay(experiments(oc))$SI[2,])))
       if (winsorizeSI) newdf$SI = ifelse(newdf$SI < 0,
                 min(newdf$SI[newdf$SI > 0]), newdf$SI)
       newdf$text = 
          as.character(paste0("ID=",newdf$id,
                          "<br>Mats=", round(newdf$mats120,3),
                          "<br>SI=", round(newdf$SI,6)))
       ggplotly(ggplot(newdf, aes(x=mats120, y=SI, 
                                colour=converged, text=text)) + 
                   geom_point() + scale_x_sqrt() + xlab("mats120") +
                   scale_y_log10(), tooltip="text")
       } )
    output$demo = renderPlot( {
       fit1 = minmodByID(oc, input$idpick, ...)
       print(plot_OGTT_fit(fit1))
       } )
    output$biplots_glucose = renderPlotly( {
       CH1 = gsub("%%N%%", input$pcbipick_1, "PC%%N%%")
       CH2 = gsub("%%N%%", input$pcbipick_2, "PC%%N%%")
       choices = as.numeric(c(input$pcbipick_1, input$pcbipick_2))
       ggplotly(ggbiplot(pglu, choices = choices, labels=gluids) + xlab(CH1) + ylab(CH2) + 
           theme_gray())
       } )
    output$biplots_insulin = renderPlotly( {
       CH1 = gsub("%%N%%", input$pcbipick_1, "PC%%N%%")
       CH2 = gsub("%%N%%", input$pcbipick_2, "PC%%N%%")
       choices = as.numeric(c(input$pcbipick_1, input$pcbipick_2))
       ggplotly(ggbiplot(pins, choices = choices, labels=insids) + xlab(CH1) + ylab(CH2) + 
           theme_gray())
       } )

    }
 shinyApp(ui, server)
}
     
