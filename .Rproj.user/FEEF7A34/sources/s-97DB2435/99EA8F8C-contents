diagnoseUI <- function(id){
  # for internal namespace structure
  ns <- NS(id)
  
  # encapsulate everything in taglist, see https://shiny.rstudio.com/articles/modules.html
  tagList(
    # if(SSO@stan_used == TRUE & SSO@stan_method == "sampling" & SSO@stan_algorithm == "NUTS") 
    uiOutput(ns("HMC"))#,
    # if(SSO@stan_used == TRUE & SSO@stan_method == "sampling" & SSO@stan_algorithm != "NUTS") uiOutput(ns("MCMC")),
    # if(SSO@stan_used == TRUE & SSO@stan_method == "variational") uiOutput(ns("VI")),
    # if(SSO@stan_used == FALSE) uiOutput(ns("MCMC"))
  )
}

diagnose <- function(input, output, session){
  
    getParcoordPlot <- callModule(parallelCoordinates, "parallelCoordinates")
    getPairsPlot <- callModule(pairs, "pairs")
    getDivergentTransitionsPlot <- callModule(divergentTransitions, "divergentTransitions")
    getDivergentScatterPlot <- callModule(divergentScatter, "divergentScatter")
    getEnergyPlot <- callModule(energy, "energy")
    getTreedepthPlot <- callModule(treedepth, "treedepth")
    getStepSizePlot <- callModule(stepSize, "stepSize")
    getAcceptancePlot <- callModule(acceptance, "acceptance")
    # 
    getTracePlot <- callModule(tracePlot, "tracePlot")
    getRankPlot <- callModule(rankPlot, "rankPlot")
    getRhatNeffSEmeanPlots <- callModule(rhat_n_eff_se_mean, "rhat_n_eff_se_mean")
    getAutoCorrelationPlot <- callModule(autoCorrelation, "autoCorrelation")
    # 
    # callModule(statsTableHMC, "statsTableHMC")
    # callModule(rhat_n_eff_se_mean_stats, "rhat_n_eff_se_mean_stats")
    # 
    # getDiagnosePlots <- reactive({
    #   list("divergentScatterPlot" = getDivergentScatterPlot(),
    #        "pairsPlot" = getPairsPlot(),
    #        "parcoordPlot" = getParcoordPlot(),
    #        "divergentTransitionsPlot" = getDivergentTransitionsPlot(),
    #        "energyPlot" = getEnergyPlot(),
    #        "treedepthPlot" = getTreedepthPlot(),
    #        "stepSizePlot" = getStepSizePlot(),
    #        "acceptancePlot" = getAcceptancePlot(),
    #        "autoCorrelationPlot" = getAutoCorrelationPlot(),
    #        "tracePlot" = getTracePlot(),
    #        "rankPlot" = getRankPlot(),
    #        "rhatPlot" = getRhatNeffSEmeanPlots()["rhatPlot"],
    #        "n_effPlot" = getRhatNeffSEmeanPlots()["n_effPlot"],
    #        "se_meanPlot" = getRhatNeffSEmeanPlots()["se_meanPlot"])
    # })
    # 
    # 
    # callModule(report, "report", ggplotsList = getDiagnosePlots, reportType = "diagnose")
  
  output$HMC <- renderUI({
    tagList(
      tags$head(
        tags$script("src"="func.js")),
      tabsetPanel(
        tabPanel(
          title = "Plots",
          id = session$ns("visualHMC"),
          navlistPanel(
            id = session$ns("HMC_navlist_vis"),
            "NUTS/HMC",
            tabPanel(
              title = "Divergent Scatter",
              id = session$ns("divergentScatterTab"),
              divergentScatterUI(session$ns("divergentScatter"))
            ),
            tabPanel(
              title = "Parallel Coordinates",
              id = session$ns("parallelCoordinatesTab"),
              parallelCoordinatesUI(session$ns("parallelCoordinates"))
            ),
            tabPanel(
              title = "Pairs",
              id = session$ns("pairsTab"),
              pairsUI(session$ns("pairs"))
            ),
            tabPanel(
              title = "Divergence Information",
              id = session$ns("divergentTransitionsTab"),
              divergentTransitionsUI(session$ns("divergentTransitions"))
            ),
            tabPanel(
              title = "Energy Information",
              id = session$ns("energyTab"),
              energyUI(session$ns("energy"))
            ),
            tabPanel(
              title = "Treedepth Information",
              id = session$ns("treedepthTab"),
              treedepthUI(session$ns("treedepth"))
            ),
            tabPanel(
              title = "Step Size Information",
              id = session$ns("stepSizeTab"),
              stepSizeUI(session$ns("stepSize"))
            ),
            tabPanel(
              title = "Acceptance Information",
              id = session$ns("acceptanceTab"),
              acceptanceUI(session$ns("acceptance"))
            ),
            "MCMC",
            tabPanel(
              title = "Autocorrelation",
              id = session$ns("autocorrelationTab"),
              autoCorrelationUI(session$ns("autoCorrelation"))
            ),
            tabPanel(
              title = "Rank Plots",
              id = session$ns("rankTab"),
              rankPlotUI(session$ns("rankPlot"))
            ),
            tabPanel(
              title = withMathJax("\\(\\hat{R}, \\text{ } n_{eff}, \\text{ se}_{mean}\\)"),
              id = session$ns("rhat_n_eff_se_meanTab"),
              value = "rhat_neff_se_mean_plot_tab",
              rhat_n_eff_se_meanUI(session$ns("rhat_n_eff_se_mean"))
            ),
            tabPanel(
              title = "Trace Plots",
              id = session$ns("traceTab"),
              tracePlotUI(session$ns("tracePlot"))
            )
          )
        ),
        tabPanel(
          title = "Stats",
          id = session$ns("numericalHMC"),
          navlistPanel(
            id = session$ns("HMC_navlist_num"),
            "NUTS/HMC")
          )#,
        # tabPanel(
        #   title = "Report",
        #   reportUI(session$ns("report"))
        # )
      )
    )
  })
  
 
}
