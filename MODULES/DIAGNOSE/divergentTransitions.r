divergentTransitionsUI <- function(id){
  ns <- NS(id)
  
  tagList(
    wellPanel(
      fluidRow(
        column(width = 4), 
        column(width = 4),
        column(width = 4, align = "right",
                 div(style = "width: 100px;",
                     numericInput(
                       ns("diagnostic_chain"),
                       label = h5(textOutput(ns("diagnostic_chain_text"))),
                       value = 0,
                       min = 0,
                       # don't allow changing chains if only 1 chain
                       max = ifelse(sso@n_chain == 1, 0, sso@n_chain)
                     )
                 )
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


divergentTransitions <- function(input, output, session){
  
  visualOptions <- callModule(plotOptions, "options")
  
  
    chain <- reactive(input$diagnostic_chain)
    include <- reactive(input$report)
    
    observe({
      toggle("caption", condition = input$showCaption)
    })
    
    output$diagnostic_chain_text <- renderText({
      if (chain() == 0)
        return("All chains")
      paste("Chain", chain())
    })
    
  
  plotOut <- function(chain){
    
    if(chain != 0) {
      mcmc_nuts_divergence(
        x = nuts_params(list(sso@sampler_params[[chain]]) %>%
                          lapply(., as.data.frame) %>%
                          lapply(., filter, row_number() > sso@n_warmup) %>%
                          lapply(., as.matrix)),
        lp = data.frame(Iteration = rep(1:(sso@n_iter - sso@n_warmup), 1),
                        Value = c(sso@posterior_sample[(sso@n_warmup + 1):sso@n_iter, chain,"log-posterior"]),
                        Chain = rep(chain, each = (sso@n_iter - sso@n_warmup))) 
      )
    } else {
      mcmc_nuts_divergence(
        x = nuts_params(sso@sampler_params %>%
                          lapply(., as.data.frame) %>%
                          lapply(., filter, row_number() > sso@n_warmup) %>%
                          lapply(., as.matrix)),
        lp = data.frame(Iteration = rep(1:(sso@n_iter - sso@n_warmup), sso@n_chain),
                        Value = c(sso@posterior_sample[(sso@n_warmup + 1):sso@n_iter, ,"log-posterior"]),
                        Chain = rep(1:sso@n_chain, each = (sso@n_iter - sso@n_warmup))) 
      )
    }
  }
  
  
  output$plot1 <- renderPlot({
    save_old_theme <- bayesplot_theme_get()
    color_scheme_set(visualOptions()$color)
    bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
    out <- plotOut(chain = chain()) 
    bayesplot_theme_set(save_old_theme)
    out
    
  })
  
  captionOut <- function(){
    HTML(paste0("These are plots of the <i> divergent transition status</i> (x-axis)",
                " against the <i> log-posterior </i> (y-axis top panel) and against the",
                " <i> acceptance statistic </i> (y-axis bottom panel) of the sampling algorithm for ",
                tolower(if (chain() == 0) {"All chains"} else {paste("Chain", chain())}), ".",
                " Divergent transitions can indicate problems for the validity of the results.",
                " A good plot would show no divergent transitions.",
                " If the divergent transitions show the same pattern as the non divergent transitions,",
                " this could indicate that the divergent transitions are false positives.",
                " A bad plot would shows systematic differences between the divergent transitions and",
                " non-divergent transitions.",
                " For more information see ",
                tags$a('https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup'), "."))
  }
  output$caption <- renderUI({
    captionOut()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'divergentTransitionPlot.pdf',
    content = function(file) {
      # ggsave(file, gridExtra::arrangeGrob(grobs = downloadSelection()))
      pdf(file)
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- plotOut(chain = chain()) 
      bayesplot_theme_set(save_old_theme)
      print(out)
      dev.off()
    })
  
  
  output$downloadRDS <- downloadHandler(
    filename = 'divergentTransitionPlot.rds',
    content = function(file) {
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- plotOut(chain = chain()) 
      bayesplot_theme_set(save_old_theme)
      saveRDS(out, file)
    })  
  
  return(reactive({
    if(include() == TRUE){
      # customized plot options return without setting the options for the other plots
      save_old_theme <- bayesplot_theme_get()
      color_scheme_set(visualOptions()$color)
      bayesplot_theme_set(eval(parse(text = select_theme(visualOptions()$theme)))) 
      out <- list(plot = plotOut(chain = chain()), 
                  caption = captionOut())
      bayesplot_theme_set(save_old_theme)
      out
    } else {
      NULL
    }
  }))
  
}