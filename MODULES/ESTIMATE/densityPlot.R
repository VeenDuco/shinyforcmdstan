densityPlotUI <- function(id){
  ns <- NS(id)
  tagList(
    wellPanel(
      fluidRow(
        column(width = 6, 
               verticalLayout(
                 selectizeInput(
                   inputId = ns("diagnostic_param"),
                   label = h5("Parameters"),
                   multiple = TRUE,
                   choices = .make_param_list_with_groups(sso),
                   selected = c(sso@param_names[1])
                 )
               )
        ),
        column(width = 4),
        column(width = 2, align = "right"
        )
               ),
      fluidRow(
        align = "right",
        plotOptionsUI(ns("options"))
      )
      ),
    plotOutput(ns("plot1")),
    checkboxInput(ns("showCaption"), "Show Caption", value = TRUE),
    hidden(
      uiOutput(ns("caption"))
    ),
    hr(), 
    # checkboxInput(ns("report"), "Include in report?")
    downloadButton(ns('downloadPlot'), 'Download Plot', class = "downloadReport"),
    downloadButton(ns('downloadRDS'), 'Download RDS', class = "downloadReport")
  )
}



densityPlot <- function(input, output, session){
  
  visualOptions <- callModule(plotOptions, "options")
  
  param <- debounce(reactive(unique(.update_params_with_groups(params = input$diagnostic_param,
                                                               all_param_names = sso@param_names))),
                    500)
  include <- reactive(input$report)
  
  observe({
    toggle("caption", condition = input$showCaption)
  })
  
  plotOut <- function(parameters, chain){
    
    validate(
      need(length(param()) > 0, "Select at least one parameter.")
    )
    mcmc_dens(
      sso@posterior_sample[(1 + sso@n_warmup) : sso@n_iter, , ],
      pars = parameters
    )
  }
  
  output$plot1 <- renderPlot({
    save_old_theme <- bayesplot_theme_get()
    color_scheme_set(visualOptions()$color)
    bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
    out <- plotOut(parameters = param())
    bayesplot_theme_set(save_old_theme)
    out
  })
  
  captionOut <- function(parameters){
    # HTML(paste0(if(length(parameters) == 1) {"This is a density plot of <i>"} else {"These are density plots of <i>"}, 
    #             paste(parameters[1:(length(parameters)-1)], collapse = ", "),
    #             if(length(parameters) > 1) {"</i> and <i>"}, 
    #             if(length(parameters) > 1) {parameters[length(parameters)]},"</i>", "."
    #             ))
    HTML(paste0(if(length(parameters) == 1) {"This is a posterior density plot."} else {"These are posterior density plots."}))
  }
  output$caption <- renderUI({
    captionOut(parameters = param())
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'densityPlot.pdf',
    content = function(file) {
      # ggsave(file, gridExtra::arrangeGrob(grobs = downloadSelection()))
      pdf(file)
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- plotOut(parameters = param())
      bayesplot_theme_set(save_old_theme)
      print(out)
      dev.off()
    })
  
  
  output$downloadRDS <- downloadHandler(
    filename = 'densityPlot.rds',
    content = function(file) {
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- plotOut(parameters = param())
      bayesplot_theme_set(save_old_theme)
      saveRDS(out, file)
    })  
  
  return(reactive({
    if(include() == TRUE){
      # customized plot options return without setting the options for the other plots
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- list(plot = plotOut(parameters = param()),
                  caption = captionOut(parameters = param()))
      bayesplot_theme_set(save_old_theme)
      out
    } else {
      NULL
    }
  }))
  
}