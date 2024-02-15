library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(svglite)

# Define UI
ui <- fluidPage(
  titlePanel("ANCOMBC Visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      uiOutput("lfc_dropdown"),
      uiOutput("se_dropdown"),
      sliderInput("threshold", "Threshold for LFC", 
                  min = 0, max = 2, value = 0.25, step = 0.05),
      actionButton("plot_btn", "Plot"),
      downloadButton("downloadPlot", "Download Plot"),
      verbatimTextOutput("info_text")
    ),
    mainPanel(
      plotOutput("waterfall_plot")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Load and prepare data
  data <- reactive({
    req(input$file1)
    df <- read_csv(input$file1$datapath)
    return(df)
  })
  
  output$lfc_dropdown <- renderUI({
    df <- data()
    selectInput("lfc_column", "Select LFC Column", choices = colnames(df))
  })
  
  output$se_dropdown <- renderUI({
    df <- data()
    selectInput("se_column", "Select Standard Error Column", choices = colnames(df))
  })
  
  output$waterfall_plot <- renderPlot({
    req(input$lfc_column, input$se_column)
    df <- data()
    lfc_column <- sym(input$lfc_column)
    se_column <- sym(input$se_column)
    
    df_filtered <- df %>%
      filter(!!lfc_column >= input$threshold | !!lfc_column <= -input$threshold) %>%
      arrange(desc(!!lfc_column)) %>%
      mutate(direct = ifelse(!!lfc_column > 0, "Positive LFC", "Negative LFC"),
             color = ifelse(!!lfc_column >= 0.25 | !!lfc_column <= -0.25, "aquamarine3", "black"))
    
    ggplot(df_filtered, aes(x = reorder(taxon, !!lfc_column), y = !!lfc_column, fill = direct)) +
      geom_col(color = "black") +
      geom_errorbar(aes(ymin = !!lfc_column - !!se_column, ymax = !!lfc_column + !!se_column), width = 0.4) +
      scale_fill_manual(values = c("Positive LFC" = "blue", "Negative LFC" = "red")) +
      labs(x = NULL, y = "Log Fold Change (LFC)", title = "Waterfall Plot of LFC") +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "top")
  })
  
  # Display information text
  output$info_text <- renderPrint({
    "Blue color represents taxa that are differentially abundant in the treated group."
  })
  
  # Download functionality
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("waterfall_plot", Sys.Date(), ".eps", sep = "")
    },
    content = function(file) {
      g <- ggplotGrob(plotOutput()$waterfall_plot())
      postscript(file, width = 10, height = 8)
      grid::grid.draw(g)
      dev.off()
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
