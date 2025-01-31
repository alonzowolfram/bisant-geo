#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(NanoStringNCTools) # For NanoString stuff.
library(GeomxTools) # For NanoString GeoMx stuff. 

# Define UI for application that draws a table.
ui <- fluidPage(
  
  # Application title
  titlePanel("QC: Probes"),
  
  # Sidebar for segment gene detection.
  sidebarLayout(
    sidebarPanel(
      # Minimum number of genes detected (minGenesDetected, min_genes_detected)
      numericInput("min_genes_detected",
                   "Minimum number of genes detected",
                   min = 0,
                   max = 20000,
                   value = 400),
      # LOQ settings
      numericInput("cutoff", "Limit of quantification (LOQ) cutoff", value = 2, min = 0, max = 100), # Limit of quantification (LOQ) cutoff, i.e. the number of SDs above mean negative probe expression a probe must be in order to be included
      numericInput("min_loq", "Minimum limit of quantification", value = 2, min = 0, max = 100) # Minimum limit of quantification (LOQ) to which values falling below cutoffLOQ will be thresholded.
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("static")
    )
  ),
  
  # Sidebar for probe gene detection.
  sidebarLayout(
    sidebarPanel(
      # Minimum number of genes detected (minGenesDetected, min_genes_detected)
      numericInput("min_gene_detection_rate",
                   "Minimum detection rate (%) for a gene to be kept",
                   min = 0,
                   max = 100,
                   value = 5)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("static2")
    )
  )
)

# Define server logic required to make a table.
server <- function(input, output) {
  
  output$static <- renderTable({
    # Load in the data object.
    data_object <- readRDS("Rdata/NanoStringGeoMxSet_qc-segments_main-module.rds")
    
    # Access the PKC files, to ensure that expected PKCs have been loaded for this study.
    pkcs <- annotation(data_object)
    modules <- base::gsub(".pkc", "", pkcs)
    pkc_summary <- data.frame(PKCs = pkcs, modules = modules)
    
    # Collapse to targets.
    target_data_object <- aggregateCounts(data_object)
    
    # Calculate LOQ per module tested.
    loq <- data.frame(row.names = colnames(target_data_object))
    for(module in modules) {
      vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                     module)
      if(all(vars[1:2] %in% colnames(pData(target_data_object)))) {
        loq[, module] <-
          pmax(input$min_loq,
               pData(target_data_object)[, vars[1]] * 
                 pData(target_data_object)[, vars[2]] ^ input$cutoff)
      }
    }
    pData(target_data_object)$LOQ <- loq
    
    # Filter out either segments and/or genes with abnormally low signal. 
    loq_mat <- c()
    for(module in modules) {
      ind <- fData(target_data_object)$Module == module
      Mat_i <- t(esApply(target_data_object[ind, ], MARGIN = 1,
                         FUN = function(x) {
                           x > loq[, module]
                         }))
      loq_mat <- rbind(loq_mat, Mat_i)
    }
    # Ensure ordering since this is stored outside of the geomxSet.
    loq_mat <- loq_mat[fData(target_data_object)$TargetName, ]
    
    # Segment gene detection:
    # Calculate detection rate information and save to phenodata.
    pData(target_data_object)$GenesDetected <- 
      colSums(loq_mat, na.rm = TRUE)
    pData(target_data_object)$GeneDetectionRate <-
      pData(target_data_object)$GenesDetected / nrow(target_data_object)
    
    # Gene detection rate:
    # Calculate detection rate:
    loq_mat <- loq_mat[, colnames(target_data_object)]
    fData(target_data_object)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
    fData(target_data_object)$DetectionRate <-
      fData(target_data_object)$DetectedSegments / nrow(pData(target_data_object))
    
    # Generate the QC summary tables. 
    # Segment gene detection.
    n_segments <- sum(pData(target_data_object)$GenesDetected >= input$min_genes_detected)
    corresponding_gene_detection_rate <- paste0(round(100*(input$min_genes_detected/nrow(target_data_object))), "%")
    qc_summary <- data.frame(
      `Number of segments (Pass)` = n_segments,
      `Number of segments (Fail)` = nrow(pData(target_data_object)) - n_segments,
      `Gene detection rate per segment` = corresponding_gene_detection_rate
    )
    # Gene detection rate 
    
    # Display the QC summary table.
    qc_summary
    
  })
  
  output$static2 <- renderTable({
    # Load in the data object.
    data_object <- readRDS("Rdata/NanoStringGeoMxSet_qc-segments_main-module.rds")
    
    # Access the PKC files, to ensure that expected PKCs have been loaded for this study.
    pkcs <- annotation(data_object)
    modules <- base::gsub(".pkc", "", pkcs)
    pkc_summary <- data.frame(PKCs = pkcs, modules = modules)
    
    # Collapse to targets.
    target_data_object <- aggregateCounts(data_object)
    
    # Calculate LOQ per module tested.
    loq <- data.frame(row.names = colnames(target_data_object))
    for(module in modules) {
      vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                     module)
      if(all(vars[1:2] %in% colnames(pData(target_data_object)))) {
        loq[, module] <-
          pmax(input$min_loq,
               pData(target_data_object)[, vars[1]] * 
                 pData(target_data_object)[, vars[2]] ^ input$cutoff)
      }
    }
    pData(target_data_object)$LOQ <- loq
    
    # Filter out either segments and/or genes with abnormally low signal. 
    loq_mat <- c()
    for(module in modules) {
      ind <- fData(target_data_object)$Module == module
      Mat_i <- t(esApply(target_data_object[ind, ], MARGIN = 1,
                         FUN = function(x) {
                           x > loq[, module]
                         }))
      loq_mat <- rbind(loq_mat, Mat_i)
    }
    # Ensure ordering since this is stored outside of the geomxSet.
    loq_mat <- loq_mat[fData(target_data_object)$TargetName, ]
    
    # Gene detection rate:
    # Calculate detection rate:
    loq_mat <- loq_mat[, colnames(target_data_object)]
    fData(target_data_object)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
    fData(target_data_object)$DetectionRate <-
      fData(target_data_object)$DetectedSegments / nrow(pData(target_data_object))
    
    # Generate the QC summary tables. 
    n_genes <- sum(fData(target_data_object)$DetectionRate >= (input$min_gene_detection_rate)/100) # Number of genes that pass detection rate cutoff. 
    qc_summary_gene <- data.frame(
      `Number of genes (Pass)` = n_genes,
      `Number of genes (Fail)` = nrow(fData(target_data_object)) - n_genes
    )
    
    # Display the QC summary table.
    qc_summary_gene
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
