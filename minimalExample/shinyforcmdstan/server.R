#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {


    
    StanFit <- reactive({
        file <- input$cmdstanInput
        ext <- tools::file_ext(file$datapath)

        req(file)
        validate(need(ext == "csv", "Please upload a csv file"))

        stanfit <- rstan::read_stan_csv(file$datapath)
        stanfit
        }
        )
    
    Draws <- reactive({
        posterior::as_draws(StanFit())
    })
    
    output$plot <- renderPlot({
        bayesplot::mcmc_dens(Draws())
    })
    

    output$modelname <- renderText({
        StanFit()@model_name
        })
    
    output$table <- renderTable({
        tt <- posterior::as_draws(StanFit())
        tt
    })
    
    output$contents <- renderTable({
       StanFit()
    })


})
