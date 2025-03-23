library(shiny)
library(ggplot2)
library(gridExtra)
library(lpirfs)
library(grid)

source("utils.r")
source("load_data.r")

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Local Projection"),
  sidebarLayout(
    sidebarPanel(
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
    # Sample data
    data <- df_ERPT_Arg  # Replace with your dataset
    endog <- c("E", "ipc", "pbird")
    exog <- c("impp_usa")
    
    results <- run_lp_model(
      data, endog, exog, max_lags = 2, 
      signif = input$signif, cumulative = input$cumulative
    )
    
    plots <- plot_lin(results)  # List of ggplot objects
    
    # Extract titles from ggplot objects
    plot_titles <- lapply(plots, function(p) {
      p$labels$title  # Extract title text
    })

    num_cols <- length(endog)  # Number of columns
    num_rows <- ceiling(length(plots) / num_cols)  # Number of rows
    
    # Arrange plots with titles
    grid.arrange(
      grobs = mapply(function(p, title) arrangeGrob(p, top = textGrob(title, gp = gpar(fontsize = 12, fontface = "bold"))), 
                     plots, plot_titles, SIMPLIFY = FALSE),
      ncol = num_cols, 
      nrow = num_rows,
      top = textGrob("Local Projection Results", gp = gpar(fontsize = 16, fontface = "bold"))
    )
  })
}


# Run the Shiny App
shinyApp(ui = ui, server = server)
