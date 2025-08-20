library(shiny)
library(ggplot2)
library(gridExtra)
library(grid)
library(lpirfs)

source("../utils.r")
source("load_data.r")

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Local Projection"),
  sidebarLayout(
    sidebarPanel(
      textInput("endog", "Endogenous Variables (comma-separated)", "E, ipc, pbird"),
      textInput("exog", "Exogenous Variables (comma-separated)", "impp_usa"),
      numericInput("max_lags", "Max Lags", value = 3, min = 1),
      textInput("newey_lags", "Newey Lags", ""),  # Allow empty input
      numericInput("horizons", "Horizons", value = 10, min = 1),
      selectInput("signif", "Significance Level", choices = c(0.95, 0.68), selected = 0.95),
      checkboxInput("cumulative", "Cumulative", FALSE)
    ),
    mainPanel(
      plotOutput("lp_plot", height = "800px")  # Increased height for better grid layout
    )
  )
)

server <- function(input, output) {
  output$lp_plot <- renderPlot({
    # Convert text inputs to character vectors, or set to NULL if empty
    endog <- if (nchar(input$endog) > 0) strsplit(input$endog, ",\\s*")[[1]] else NULL
    exog <- if (nchar(input$exog) > 0) strsplit(input$exog, ",\\s*")[[1]] else NULL
    newey_lags <- suppressWarnings(as.integer(input$newey_lags))
    if (is.na(newey_lags)) newey_lags <- NULL

    # Sample data
    data <- df_ERPT_Arg  # Replace with your dataset
    
    results <- run_lp_model(
      data, 
      endog, 
      exog, 
      max_lags = input$max_lags, 
      newey_lags = newey_lags, 
      horizons = input$horizons,
      signif = input$signif, 
      cumulative = input$cumulative
    )
    
    plots <- plot_lin(results)  # List of ggplot objects
    
    # Extract titles from ggplot objects
    plot_titles <- lapply(plots, function(p) {
      p$labels$title  # Extract title text
    })

    num_cols <- length(endog)  # Number of columns
    num_rows <- ceiling(length(plots) / num_cols)  # Number of rows
    
    top_title <- if (input$cumulative) {
      "Cumulative Local Projection Results"
    } else {
      "Local Projection Results"
    }

    # Arrange plots with titles
    grid.arrange(
      grobs = mapply(function(p, title) arrangeGrob(p, top = textGrob(title, gp = gpar(fontsize = 12, fontface = "bold"))), 
                     plots, plot_titles, SIMPLIFY = FALSE),
      ncol = num_cols, 
      nrow = num_rows,
      top = textGrob(top_title, gp = gpar(fontsize = 16, fontface = "bold"))
    )
  })
}

# Run the Shiny App
shinyApp(ui = ui, server = server)