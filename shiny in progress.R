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
                 fileInput("metadata", "Metadata", accept = c(".csv", ".tsv", ".txt")),
                 fileInput("expression_matrix", "Expression Matrix", accept = c(".csv", ".tsv", ".txt")),
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
                 h3("Annotation"),
                 fileInput("annotation_model", "Upload Annotation Model", accept = c(".csv", ".tsv", ".txt")),
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
  rv <- reactiveValues(metadata = NULL, expression_matrix = NULL, seurat_object = NULL, annotation_model = NULL, dge_results = NULL)
  
  # Load Data when "Load Data" button is clicked
  observeEvent(input$load_data, {
    # Load Metadata
    if (!is.null(input$metadata)) {
      rv$metadata <- read.csv(input$metadata$datapath)
    }
    
    # Load Expression Matrix
    if (!is.null(input$expression_matrix)) {
      rv$expression_matrix <- read.csv(input$expression_matrix$datapath, row.names = 1)
    }
    
    # Load Seurat Object
    if (!is.null(input$seurat_object)) {
      rv$seurat_object <- readRDS(input$seurat_object$datapath)
      
      # Check for meta.data in the Seurat object
      if (is.null(rv$seurat_object@meta.data) || nrow(rv$seurat_object@meta.data) == 0) {
        showNotification("Seurat object does not contain metadata!", type = "error")
      }
    }
    
    # Update Previews
    output$metadata_preview <- renderDT({
      req(rv$metadata)
      datatable(rv$metadata, options = list(pageLength = 5))
    })
    
    output$expression_matrix_preview <- renderDT({
      req(rv$expression_matrix)
      datatable(rv$expression_matrix, options = list(pageLength = 5))
    })
    
    output$seurat_metadata_preview <- renderDT({
      req(rv$seurat_object)
      req(rv$seurat_object@meta.data)
      datatable(rv$seurat_object@meta.data, options = list(pageLength = 5))
    })
    
    showNotification("Data loaded successfully!", type = "message")
  })
  
  # Reset Button Logic
  observeEvent(input$reset, {
    rv$metadata <- NULL
    rv$expression_matrix <- NULL
    rv$seurat_object <- NULL
    rv$annotation_model <- NULL
    rv$dge_results <- NULL
    output$metadata_preview <- renderDT(NULL)
    output$expression_matrix_preview <- renderDT(NULL)
    output$seurat_metadata_preview <- renderDT(NULL)
    showNotification("All data has been reset.", type = "warning")
  })
  
  # Perform Annotation
  observeEvent(input$perform_annotation, {
    req(rv$seurat_object)
    
    if (!is.null(input$annotation_model)) {
      file_type <- file_ext(input$annotation_model$name)
      
      if (file_type == "rds") {
        rv$annotation_model <- readRDS(input$annotation_model$datapath)
        rv$seurat_object@meta.data$predicted_annotations <- predict(rv$annotation_model, newdata = t(rv$seurat_object@assays$RNA@data))
      } else if (file_type == "csv") {
        rules <- read.csv(input$annotation_model$datapath)
        rv$seurat_object@meta.data$predicted_annotations <- apply(
          rv$seurat_object@meta.data, 1, function(row) {
            match <- rules$label[rules$marker %in% row]
            ifelse(length(match) > 0, match[1], "Unlabeled")
          }
        )
      }
    }
    
    output$annotation_table <- renderDT({
      req(rv$seurat_object@meta.data)
      datatable(rv$seurat_object@meta.data, options = list(pageLength = 5))
    })
    
    output$annotation_barplot <- renderPlot({
      req(rv$seurat_object@meta.data$predicted_annotations)
      ggplot(as.data.frame(rv$seurat_object@meta.data), aes(x = predicted_annotations)) +
        geom_bar(fill = "skyblue") +
        theme_minimal() +
        labs(title = "Cell Type Distribution", x = "Cell Type", y = "Count")
    })
  })
  
  # Run DGE Analysis
  observeEvent(input$run_dge, {
    req(input$group_column, input$group1, input$group2, rv$seurat_object)
    
    group1_cells <- which(rv$seurat_object@meta.data[[input$group_column]] == input$group1)
    group2_cells <- which(rv$seurat_object@meta.data[[input$group_column]] == input$group2)
    
    rv$dge_results <- FindMarkers(rv$seurat_object, ident.1 = group1_cells, ident.2 = group2_cells)
    
    output$dge_table <- renderDT({
      req(rv$dge_results)
      datatable(rv$dge_results, options = list(pageLength = 5))
    })
    
    output$volcano_plot <- renderPlot({
      req(rv$dge_results)
      ggplot(rv$dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = avg_log2FC > 0), size = 1) +
        theme_minimal() +
        labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-Value")
    })
    
    output$dge_heatmap <- renderPlot({
      req(rv$dge_results)
      top_genes <- rownames(head(rv$dge_results[order(rv$dge_results$p_val_adj), ], 20))
      heatmap_data <- as.matrix(rv$seurat_object@assays$RNA@data[top_genes, ])
      pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, main = "Top DGE Heatmap")
    })
    
    output$violin_plot <- renderPlot({
      req(input$gene_of_interest)
      VlnPlot(rv$seurat_object, features = input$gene_of_interest, group.by = input$group_column)
    })
  })
}

# Create the app

shinyApp(ui = ui, server = server)
