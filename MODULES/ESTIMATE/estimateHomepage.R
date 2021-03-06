estimateUI <- function(id){
  # for internal namespace structure
  ns <- NS(id)
  
  # encapsulate everything in taglist, see https://shiny.rstudio.com/articles/modules.html
  tagList(
    uiOutput(ns("estimateHomepage"))
  )
  
}

estimate <- function(input, output, session){
  
  getVisualPlots <- callModule(visualEstimate, "visualEstimate")
  callModule(numericalEstimate, "numericalEstimate")

  # getEstimatePlots <- reactive({
  #   list("intervalsPlot" = getVisualPlots()["intervalsPlot"],
  #        "areasPlot" = getVisualPlots()["areasPlot"],
  #        "scatterPlot" = getVisualPlots()["scatterPlot"],
  #        "densityPlot" = getVisualPlots()["densityPlot"],
  #        "histogramPlot" = getVisualPlots()["histogramPlot"]
  #   )
  # })
  # 
  # callModule(report, "report", ggplotsList = getEstimatePlots, reportType = "estimate")
  
  output$estimateHomepage <- renderUI({
    tagList(
      tabsetPanel(
        id = session$ns("diagnose_tabset"),
        tabPanel(
          title = "Plots",
          id = session$ns("visualEstimateTab"),
        visualEstimateUI(session$ns("visualEstimate"))
        ),
        tabPanel(
          title = "Stats",
          id = session$ns("numericalEstimateTab"),
        numericalEstimateUI(session$ns("numericalEstimate"))
        )#,
        # tabPanel(
        #   title = "Report",
        #   id = session$ns("reportTab"),
        #   reportUI(session$ns("report"))
        # )
      )
    )
  })
}
