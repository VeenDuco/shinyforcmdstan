#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            fileInput("cmdstanInput", "Upload your csv file from cmdstan",
                      multiple = FALSE,
                      accept = ".csv"),
            # actionButton("processInput", "Process file"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            textOutput("modelname"),
            plotOutput("plot"),
            tableOutput("table"),
            tableOutput("contents")
        )
    )
))
