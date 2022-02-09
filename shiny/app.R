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
library(DiagrammeR)
library(shinybusy)
#library(knitr)

source("scripts/gen_starts.R")
source("scripts/clpm10_c.R")
source("scripts/ri_clpm10_c.R")
source("scripts/starts_c.R")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Generating and Modeling CLPM/RI-CLPM/STARTS Data"),
    use_busy_spinner(spin="fading-circle"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("n",
                        "Total N:",
                        min = 0,
                        value = 500,
                        width = "25%"),
            h3("Variances (Using STARTS Terminology)"),
            h4("X Variable"),
            fluidRow(
                column(4, numericInput("st_x",
                                       "Stable Trait",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("ar_x",
                                       "Autoregressive",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("state_x",
                                       "State",
                                       min = 0,
                                       step = .05,
                                       value = 1))),
            h4("Y Variable"),
            fluidRow(
                column(4, numericInput("st_y",
                                       "Stable Trait",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("ar_y",
                                       "Autoregressive",
                                       min = 0,
                                       step = .05,
                                       value = 1)),
                column(4, numericInput("state_y",
                                       "State",
                                       min = 0,
                                       step = .05,
                                       value = 1))),
            h3("Autoregressive Parameters"),
            fluidRow(
                column(4, numericInput("stability_x",
                                       "Stability of X",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .5)),
                column(4, numericInput("stability_y",
                                       "Stability of Y",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .5))),
            fluidRow(
                column(4, numericInput("XonY",
                                       "Cross-lag: Y Predicting X",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .2)),
                column(4, numericInput("YonX",
                                       "Cross-lag: X Predicting Y",
                                       min = 0,
                                       max = 1,
                                       step = .05,
                                       value = .2))),
;            h3("Correlations"),
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
            actionButton("update", "Simulate Data and Run Models"),
            p(),
            h3("How To Use This App"),
            p("This Shiny App generates 10 waves of data for two variables, based on a bivariate \"Stable Trait, Autoregressive Trait, State\" (STARTS) Model. The STARTS model decomposes the variance into a Stable Trait, an Autoregressive Trait, and a State, and the App allows you to specify the variance for each component. You can also specify correlations among components, the stability of the autoregressive part, and cross-lagged paths. A five-wave version of this model is presented in the figure."),
            p("The Random Intercept Cross-Lagged Panel Model (RI-CLPM) is equivalent to a STARTS model without the state variance; the Cross-Lagged Panel Model (CLPM) is equivalent to the RI-CLPM without the stable trait. Therefore you can specify different data generating models by setting these varianes to zero."),
            p("The app then runs the CLPM, RI-CLPM, and (optionally) STARTS models using the generated data"),
            p("The purpose of this app is to test what happens to cross-lagged paths in the CLPM or RI-CLPM when either the Stable Trait or State variances are omitted from the model."),
            p("By default, the app only runs the CLPM and RI-CLPM. You can also try the full STARTS model on the generated data, but it can take a few extra seconds to run, so you have to explicitly ask for this output. Also note that you can specify values for which the data are impossible to generate (e.g., high stability plus strong cross-lagged paths). I don't have good error checking for this yet, so if the page freezes, just reload to restart.")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          textOutput("clpm_results"),
          tableOutput("clpm_table"),
          textOutput("riclpm_results"),
          tableOutput("riclpm_table"),
          textOutput("starts_results"),
          tableOutput("starts_table"),
          grVizOutput('starts', width = "50%")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$starts <- renderGrViz({
        grViz( " digraph D {
  label = 'The Stable Trait, Autoregressive Trait, State (STARTS) Model'
  ranksep=.5;
  node [shape = ellipse];
  STX [label='X Stable\nTrait']; 
  STY [label='Y Stable\nTrait'];
  node [shape = ellipse, width = 1];
  ARX1 ARX2 ARX3 ARX4 ARX5;
  ARY1 ARY2 ARY3 ARY4 ARY5;
  node [shape = box, width = 1];
  X1 X2 X3 X4 X5
  Y1 Y2 Y3 Y4 Y5
  node [shape = circle, width = .05, fontsize='10'];
  sx1 sx2 sx3 sx4 sx5
  sy1 sy2 sy3 sy4 sy5
    
  {rank = same X1 X2 X3 X4 X5}
  {rank = same ARX1 ARX2 ARX3 ARX4 ARX5}
  {rank = same ARY1 ARY2 ARY3 ARY4 ARY5}
  {rank = same Y1 Y2 Y3 Y4 Y5}
  {rank = same sx1 sx2 sx3 sx4 sx5}
  {rank = same sy1 sy2 sy3 sy4 sy5}
  
  edge [dir = back, label = '1'];
  X1 -> ARX1
  X2 -> ARX2
  X3 -> ARX3
  X4 -> ARX4
  X5 -> ARX5
  
  edge [dir = ''];
  
  ARX1 -> ARX2 -> ARX3 -> ARX4 -> ARX5
  
  edge [label ='']
  ARX1 -> ARY2
  ARX2 -> ARY3
  ARX3 -> ARY4
  ARX4 -> ARY5
  ARY1 -> ARX2
  ARY2 -> ARX3
  ARY3 -> ARX4
  ARY4 -> ARX5
  
    edge [label ='1']
  ARY1 -> ARY2 -> ARY3 -> ARY4 -> ARY5
  
  ARY1 -> Y1
  ARY2 -> Y2
  ARY3 -> Y3
  ARY4 -> Y4
  ARY5 -> Y5
  
  STX -> {X1, X2, X3, X4, X5}
  
  edge [dir=back];
  {Y1, Y2, Y3, Y4, Y5} -> STY
  
  edge [dir='']
  
  sx1 -> X1
  sx2 -> X2
  sx3 -> X3
  sx4 -> X4
  sx5 -> X5
  
  edge [dir=back]
  Y1 -> sy1
  Y2 -> sy2
  Y3 -> sy3
  Y4 -> sy4
  Y5 -> sy5
  
  edge [style=invis]
  STX -> {sx1, sx2, sx3, sx4, sx5}
  {sy1, sy2, sy3, sy4, sy5} -> STY
  
  edge [style = '', dir=both, label = ''];
  ARX1 -> ARY1
  STX -> STY
  

}"
)})
  
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
                              "Y Predicting X",
                              "X Predicting Y",
                              "Variance of X",
                              "Variance of Y",
                              "Correlation Between X and Y")
    output$clpm_table <- renderTable(clpm_table, rownames = TRUE)
    fit_riclpm <- lavaan(ri_clpm10_c, data = data)
    riclpm_table <- parameterEstimates(fit_riclpm)[c(41,50,68,59,77,78,79,80,99,100),c(5:10)]
    rownames(riclpm_table) <- c("Stability of X",
                              "Stability of Y",
                              "Y Predicting X",
                              "X Predicting Y",
                              "Variance of X Random Intercept",
                              "Variance of Y Random Intercept",
                              "Variance of Autoregressive X",
                              "Variance of Autoregressive Y",
                              "Correlation Between X and Y Random Intercepts",
                              "Correlation Between X and Y Autoregressive")
    output$riclpm_table <- renderTable(riclpm_table, rownames = TRUE)

    
    if (input$run_starts==TRUE) {
        show_spinner()
        fit_starts <- lavaan(starts_c, data = data)
        hide_spinner()
        starts_table <- parameterEstimates(fit_starts)[c(41,50,68,59,77,78,79,80,99,109,119,120,121),c(5:10)]
        rownames(starts_table) <- c("Stability of X",
                                    "Stability of Y",
                                    "Y Predicting X",
                                    "X Predicting Y",
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
