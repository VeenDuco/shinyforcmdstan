rhat_n_eff_se_mean_statsUI <- function(id){
  # for internal namespace structure
  ns <- NS(id)
  tagList(
    wellPanel(
      fluidRow(
        column(width = 6,
               selectizeInput(
                 inputId = ns("diagnostic_param"),
                 label = h5("Parameter"),
                 multiple = TRUE,
                 choices = .make_param_list_with_groups(sso),
                 selected = if(length(sso@param_names) > 10) {
                   sso@param_names[order(sso@summary[, "n_eff"])[1:10]]
                 }  else {
                   sso@param_names[order(sso@summary[, "n_eff"])]
                 } 
               )
        ), 
        column(width = 4),
        column(width = 2, align = "right",
               div(style = "width: 100px;",
                   numericInput(
                     ns("sampler_digits"),
                     label = h5("Decimals"),
                     value = 2,
                     min = 0,
                     max = 10,
                     step = 1
                   )
               )
        )
      )
    ),
    DT::dataTableOutput(ns("sampler_summary"))
  )
  
}




rhat_n_eff_se_mean_stats <- function(input, output, session){
  
  param <- reactive(unique(.update_params_with_groups(params = input$diagnostic_param,
                                                      all_param_names = sso@param_names)))
  
  digits <- reactive(input$sampler_digits)
  
  MCMCtable <- reactive({
    out <- sso@summary[, c("Rhat", "n_eff", "se_mean", "sd")]
    out[, 2] <- out[, 2] / ((sso@n_iter - sso@n_warmup) * sso@n_chain)
    out[, 3] <- out[, 3] / out[, 4]
    out <- out[param(), 1:3]
    out <- cbind(out, as.matrix(sso_monitor_summary)[, c("n_eff", "Bulk_ESS", "Tail_ESS")][param(), ])
    if(length(param()) == 1) out <- matrix(out, nrow = 1); rownames(out) <- param()
    colnames(out) <- c("Rhat", "n_eff / N", "se_mean / sd", "n_eff", "Bulk_ESS", "Tail_ESS")
    out <- round(out, digits())
    out
    
  })
  
  
  output$sampler_summary <- DT::renderDataTable({
    validate(
      need(length(param()) > 0, "Select at least one parameter.")
    )
    DT::datatable({
      MCMCtable() 
    }, options = list(
      order = list(2, 'asc'),
      processing = TRUE,
      deferRender = TRUE,
      scrollX = TRUE,
      scrollY = "300px",
      scrollCollapse = TRUE,
      paging = FALSE,
      searching = TRUE,
      info = FALSE
    ))
  })
  
  
}
