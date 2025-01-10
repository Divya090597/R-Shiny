# R-Shiny 

## Structure of Shiny app ##

3 components 
             
             - User Interface (Ui.R)

             - Server function (server.R)
             
             - Shiny app function ( Fuses the UI and Server components)

The UI is the front end that accepts user input values

The server is the backend that process these input values to finally produce the output results that are finally displayed on the website.

This Shiny app is designed to handle cell functional state annotation and differential gene expression (DGE) analysis for single-cell or spatial transcriptomics data. The code includes features for uploading data, annotating cell types (models), performing DGE analysis, and downloading results.

## Features ##

Upload Metadata: Allows users to upload a file containing metadata for cells (e.g., sample IDs, grouping information etc.).

Upload Expression Matrix: Handles gene expression data (genes as rows, cells as columns).

Upload Seurat Object: Users can upload a .rds file representing a preprocessed Seurat object

*Intention*: Provide flexibility in data inputs, supporting both standard formats and advanced data objects (Seurat).


## Differential Gene Expression Analysis ##

Dynamic Group Selection: Users can select a metadata column to define groupings (e.g., treatment vs. control) and then choose two groups for comparison.

DGE Calculation: The app runs FindMarkers() from the Seurat package, comparing the two selected groups to identify differentially expressed genes.

Visualization: Results are displayed in a table and as a volcano plot (showing log fold-change vs. significance).

**Intention**: Facilitate group comparisons and provide insights into the biological significance of gene expression changes.
