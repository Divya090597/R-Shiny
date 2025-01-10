# Load the required libraries

library(shiny)
library(shinythemes)
library(Seurat)
library(DT)
library(ggplot2)
library(pheatmap)
library(tools)

# Increase file upload limit to 500MB
options(shiny.maxRequestSize = 500*1024^2)

ui <- fluidPage(
  theme = shinytheme("united"),
  navbarPage(
    "Cell Functional State Annotation & DGE Analysis",
    
    # Tab: Data Input
    tabPanel("Data Input",
             sidebarLayout(
               sidebarPanel(
                 h3("Upload Data"),
                 fileInput("metadata", "Metadata", accept = c(".csv", ".txt")),
                 fileInput("expression_matrix", "Expression Matrix", accept = c(".csv", ".txt")),
                 fileInput("seurat_object", "Seurat Object (RDS)", accept = ".rds"),
                 actionButton("load_data", "Load Data", class = "btn-primary"),
                 actionButton("reset", "Reset", class = "btn-danger")
               ),
               mainPanel(
                 h3("Data Previews"),
                 tabsetPanel(
                   tabPanel("Metadata Preview", DTOutput("metadata_preview")),
                   tabPanel("Expression Matrix Preview", DTOutput("expression_matrix_preview")),
                   tabPanel("Seurat Object Metadata", DTOutput("seurat_metadata_preview"))
                 )
               )
             )
    ),
    
    # Tab: Annotation
    tabPanel("Cell Type Annotation",
             sidebarLayout(
               sidebarPanel(
                 h3("Annotation Options"),
                 fileInput("annotation_model", "Upload Annotation Model (RDS/CSV)", accept = c(".rds", ".csv")),
                 actionButton("perform_annotation", "Annotate", class = "btn-success"),
                 downloadButton("download_annotation", "Download Annotated Data")
               ),
               mainPanel(
                 h3("Annotation Results"),
                 tabsetPanel(
                   tabPanel("Table", DTOutput("annotation_table")),
                   tabPanel("Cell Type Distribution", plotOutput("annotation_barplot"))
                 )
               )
             )
    ),
    
    # Tab: Differential Gene Expression (DGE)
    tabPanel("Differential Gene Expression",
             sidebarLayout(
               sidebarPanel(
                 h3("DGE Analysis"),
                 selectInput("group_column", "Select Group Column", choices = NULL),
                 uiOutput("group_values"),
                 textInput("gene_of_interest", "Gene for Violin Plot", value = "GeneA"),
                 actionButton("run_dge", "Run DGE Analysis", class = "btn-warning"),
                 downloadButton("download_dge", "Download DGE Results")
               ),
               mainPanel(
                 h3("DGE Visualizations"),
                 tabsetPanel(
                   tabPanel("DGE Table", DTOutput("dge_table")),
                   tabPanel("Volcano Plot", plotOutput("volcano_plot")),
                   tabPanel("Heatmap", plotOutput("dge_heatmap")),
                   tabPanel("Violin Plot", plotOutput("violin_plot"))
                 )
               )
             )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(metadata = NULL, expression_matrix = NULL, seurat_object = NULL)
  
  # Load Data when "Load Data" button is clicked
  observeEvent(input$load_data, {
    # Load Metadata
    if (!is.null(input$metadata)) {
      rv$metadata <- read.csv(input$metadata$datapath)
      print("Metadata loaded:")
      print(head(rv$metadata))  # Debugging: Check metadata content
    }
    
    # Load Expression Matrix
    if (!is.null(input$expression_matrix)) {
      rv$expression_matrix <- read.csv(input$expression_matrix$datapath, row.names = 1)
      print("Expression matrix loaded:")
      print(dim(rv$expression_matrix))  # Debugging: Check dimensions
    }
    
    # Load Seurat Object
    if (!is.null(input$seurat_object)) {
      rv$seurat_object <- readRDS(input$seurat_object$datapath)
      print("Seurat object loaded:")
      print(rv$seurat_object)  # Debugging: Check object structure
    }
    
    # Update Previews
    output$metadata_preview <- renderDT({
      req(rv$metadata)  # Ensure metadata is loaded
      datatable(rv$metadata, options = list(pageLength = 5))
    })
    
    output$expression_matrix_preview <- renderDT({
      req(rv$expression_matrix)  # Ensure expression matrix is loaded
      datatable(rv$expression_matrix, options = list(pageLength = 5))
    })
    
    output$seurat_metadata_preview <- renderDT({
      req(rv$seurat_object)  # Ensure Seurat object is loaded
      datatable(rv$seurat_object@meta.data, options = list(pageLength = 5))
    })
    
    showNotification("Data loaded successfully!", type = "message")
  })
  
  # Reset Button Logic
  observeEvent(input$reset, {
    rv$metadata <- NULL
    rv$expression_matrix <- NULL
    rv$seurat_object <- NULL
    output$metadata_preview <- renderDT(NULL)
    output$expression_matrix_preview <- renderDT(NULL)
    output$seurat_metadata_preview <- renderDT(NULL)
    showNotification("All data has been reset.", type = "warning")
  })
  
  # Placeholder for annotation and DGE logic
  # (These sections would be added as per the rest of the app logic)
}

shinyApp(ui = ui, server = server)
