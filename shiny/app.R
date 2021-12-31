#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(lavaan)
#library(knitr)

source("scripts/gen_starts.R")
source("scripts/clpm10_c.R")
source("scripts/ri_clpm10_c.R")
source("scripts/starts_c.R")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Generating CLPM/RI-CLPM/STARTS Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("n",
                        "Total N:",
                        min = 0,
                        value = 500,
                        width = "25%"),
            titlePanel("Variances (Using STARTS Terminology)"),
            fluidRow(
                column(4, numericInput("st_x",
                                       "X Stable Trait",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("ar_x",
                                       "X Autoregressive",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("state_x",
                                       "X State",
                                       min = 0,
                                       step = .05,
                                       value = 1))),
            fluidRow(
                column(4, numericInput("st_y",
                                       "Y Stable Trait",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("ar_y",
                                       "Y Autoregressive",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("state_y",
                                       "Y State",
                                       min = 0,
                                       step = .05,
                                       value = 1))),
            titlePanel("Autoregressive Parameters"),
            fluidRow(
                column(6, numericInput("stability_x",
                                       "Stability of X",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .5)),
                column(6, numericInput("stability_y",
                                       "Stability of Y",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .5))),
            fluidRow(
                column(6, numericInput("YonX",
                                       "Cross-lag: Y regressed on X:",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .2)),
                column(6, numericInput("XonY",
                                       "Cross-lag: X regressed on Y:",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .2))),
            titlePanel("Correlations"),
            fluidRow(
                column(4, numericInput("st_cor",
                                       "Correlation Between Stable Traits:",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .5)),
                column(4, numericInput("ar_cor",
                        "Correlation Between Initial AR:",
                        min = 0,
                        max = 1,
                        step = .05,
                        value = .5)),
                ),
            checkboxInput("run_starts",
                          label = "Run Starts (Takes Longer)",
                          value = FALSE),
            actionButton("update", "Update Data")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          textOutput("clpm_results"),
          tableOutput("clpm_table"),
          textOutput("riclpm_results"),
          tableOutput("riclpm_table"),
          textOutput("starts_results"),
          tableOutput("starts_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  observeEvent(input$update, {
    data <- gen_starts(n=input$n,
                       nwaves=10,   # Number of waves
                       ri_x=input$st_x,     # Random intercept variance for X
                       ri_y=input$st_y,     # Random intercept variance for Y
                       cor_i=input$st_cor,   # Correlation between intercepts (as correlation)
                       x=input$ar_x,        # AR variance for X
                       y=input$ar_y,        # AR variance for Y
                       stab_x=input$stability_x,  # Stability of X
                       stab_y=input$stability_y,  # Stability of Y
                       yx=input$YonX,      # Cross lag (Y regressed on X)
                       xy=input$XonY,      # Cross lag (X regressed on Y)
                       cor_xy=input$ar_cor,  # Correlation between X and Y (as correlation)
                       xr=input$state_x,       # Measurement error for X
                       yr=input$state_y        # Measurement error for Y
                       )
        fit_clpm <- lavaan(clpm10_c, data = data)
    clpm_table <- parameterEstimates(fit_clpm)[c(1,2,3,4,37,38,57),c(5:10)]
    rownames(clpm_table) <- c("Stability of X",
                              "Stability of Y",
                              "X Predicted by Y",
                              "Y Predicted by X",
                              "Variance of X",
                              "Variance of Y",
                              "Correlation Between X and Y")
    output$clpm_table <- renderTable(clpm_table, rownames = TRUE)
    fit_riclpm <- lavaan(ri_clpm10_c, data = data)
    riclpm_table <- parameterEstimates(fit_riclpm)[c(41,50,68,59,77,78,79,80,99,100),c(5:10)]
    rownames(riclpm_table) <- c("Stability of X",
                              "Stability of Y",
                              "X Predicted by Y",
                              "Y Predicted by X",
                              "Variance of X Random Intercept",
                              "Variance of Y Random Intercept",
                              "Variance of Autoregressive X",
                              "Variance of Autoregressive Y",
                              "Correlation Between X and Y RI",
                              "Correlation Between X and Y Autoregressive")
    output$riclpm_table <- renderTable(riclpm_table, rownames = TRUE)

    
    if (input$run_starts==TRUE) {
        fit_starts <- lavaan(starts_c, data = data)
        starts_table <- parameterEstimates(fit_starts)[c(41,50,68,59,77,78,79,80,99,109,119,120,121),c(5:10)]
        rownames(starts_table) <- c("Stability of X",
                                    "Stability of Y",
                                    "X on Y",
                                    "Y on X",
                                    "Variance of X Random Intercept",
                                    "Variance of Y Random Intercept",
                                    "Variance of Autoregressive X",
                                    "Variance of Autoregressive Y",
                                    "Variance of State X",
                                    "Variance of State Y",
                                    "Correlation Between X and Y RI",
                                    "Correlation Between X and Y Deviations",
                                    "Correlation Between X and Y Deviation Residuals")      
        output$starts_table <- renderTable(starts_table, rownames = TRUE)
        output$starts_results <- renderText({"STARTS Results"})
    } 
    output$clpm_results <- renderText({"CLPM Results"})
    output$riclpm_results <- renderText({"RI-CLPM Results"})
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
