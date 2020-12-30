statsTableHMCUI <- function(id){
  # for internal namespace structure
  ns <- NS(id)
  tagList(
    wellPanel(
      fluidRow(
        column(width = 6,
               radioButtons(
                 ns("sampler_report"),
                 label = h5("Statistic"),
                 choices = list(
                   Mean = "average",
                   SD = "sd",
                   Max = "maximum",
                   Min = "minimum"
                 ),
                 inline = TRUE
               )
               ), 
        column(width = 4,
               radioButtons(
                 ns("sampler_warmup"),
                 label = h5("Warmup"),
                 choices = list(Omit = "omit", Include = "include"),
                 inline = TRUE
               )
               ),
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




statsTableHMC <- function(input, output, session){
  
  summary_stats_sampler <- reactive({
    validate(
      need(shinystan:::.sso_env$.SHINYSTAN_OBJECT@misc$stan_algorithm %in% c("NUTS", "HMC"), message = "Only available for algorithm = NUTS or HMC"),
      need(input$sampler_warmup, message = "Loading...")
    )
    sp <- if (input$sampler_warmup == "include")
      shinystan:::.sso_env$.SHINYSTAN_OBJECT@sampler_params else shinystan:::.sso_env$.SHINYSTAN_OBJECT@sampler_params %>%
      lapply(., as.data.frame) %>%
      lapply(., filter, row_number() > shinystan:::.sso_env$.SHINYSTAN_OBJECT@n_warmup)
    
    .sampler_stuff <- function(X, param, report) {
      sapply_funs <- function(x, fun_name) {
        funs <- list(
          maxf = function(x) max(x[, param]),
          minf = function(x) min(x[, param]),
          meanf = function(x) mean(x[, param]),
          sdf = function(x) sd(x[, param])
        )
        sapply(x, FUN = funs[[fun_name]])
      }
      out <- if (report == "maximum") sapply_funs(X, "maxf") 
      else if (report == "minimum") sapply_funs(X, "minf")
      else if (report == "sd") sapply_funs(X, "sdf")
      else sapply_funs(X, "meanf")
      
      names(out) <- paste0("chain",1:length(out))
      out
    }
    
    .sampler_summary <- function(sampler_params, warmup_val,
                                 report = "average", digits = 4){ 
      
      params <- colnames(shinystan:::.sso_env$.SHINYSTAN_OBJECT@sampler_params[[1]])
      out <- sapply(params, FUN = function(p) 
        .sampler_stuff(X = sampler_params, param = p, report = report))
      
      if (length(dim(out)) > 1) { # if multiple chains
        out <- rbind("All chains" = colMeans(out), out)
        colnames(out) <- gsub("__","",colnames(out))
        out <- formatC(round(out, digits), format = 'f', digits = digits)
      } else { # if only 1 chain
        names(out) <- gsub("__.chain1", "", names(out))
        out <- round(t(out), digits)
      }
      out
    }
    
    
    do.call(
      ".sampler_summary",
      args = list(
        sampler_params  = sp,
        warmup_val      = shinystan:::.sso_env$.SHINYSTAN_OBJECT@n_warmup,
        report          = input$sampler_report,
        digits          = input$sampler_digits
      )
    )
  })
  
  output$sampler_summary <- DT::renderDataTable({
    DT::datatable({
      summary_stats_sampler()
    }, options = list(
      # rownames = FALSE,
      processing = TRUE,
      deferRender = TRUE,
      scrollX = TRUE,
      scrollY = "300px",
      scrollCollapse = TRUE,
      paging = FALSE,
      searching = FALSE,
      info = FALSE
    ))
  })
  
}
