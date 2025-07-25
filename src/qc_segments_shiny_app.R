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

# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("QC: Segments"),
  
  # Sidebar with 1) selection of the module from the data object list and 2) slider inputs for QC metric cutoffs
  sidebarLayout(
    sidebarPanel(
      # Dynamically generated selection of module from data object list
      uiOutput("dynamicModules"),
      
      # Minimum number of segment reads (minSegmentReads, min_segment_reads)
      sliderInput("min_segment_reads",
                  "Minimum number of segment reads",
                  min = 0,
                  max = 1000,
                  value = 1000),
      # Percent reads trimmed.
      sliderInput("percent_trimmed",
                  "Minimum % of reads trimmed",
                  min = 0,
                  max = 100,
                  value = 80),
      # Percent reads stitched.
      sliderInput("percent_stitched",
                  "Minimum % of reads stitched",
                  min = 0,
                  max = 100,
                  value = 80),
      # Percent reads aligned.
      sliderInput("percent_aligned",
                  "Minimum % of reads aligned",
                  min = 0,
                  max = 100,
                  value = 80),
      # Minimum sequencing saturation (%).
      sliderInput("percent_saturation",
                  "Minimum sequencing saturation (%)",
                  min = 0,
                  max = 100,
                  value = 50),
      # Minimum number of negative control counts.
      sliderInput("min_negative_count",
                  "Minimum # of negative control counts",
                  min = 0,
                  max = 10,
                  value = 5),
      # Maximum number of counts observed in No Template Control well.
      sliderInput("max_ntc_count",
                  "Maximum # of counts in No Template Control (NTC) well",
                  min = 0,
                  max = 10000,
                  value = 1000),
      # Minimum number of nuclei estimated.
      sliderInput("min_nuclei",
                  "Minimum # of nuclei estimated",
                  min = 0,
                  max = 100,
                  value = 20),
      # Minimum segment area.
      sliderInput("min_area",
                  "Minimum segment area",
                  min = 1000,
                  max = 10000,
                  value = 5000),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("static")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Load in the data object list.
  data_object_list <- readRDS("Rdata/NanoStringGeoMxSet_raw.rds")
  
  # Dynamically generate the drop-down menu with the modules of the `data_object_list`
  output$dynamicModules <- renderUI({
    selectInput(
      inputId = "module",
      label = "Select module",
      choices = as.list(names(data_object_list)),
      multiple = FALSE,
      selectize = TRUE,
      width = NULL,
      size = NULL
    )
  })
  
  # Render the QC table
  output$static <- renderTable({
    # Set the `data_object` based on the user's selection
    data_object <- data_object_list[[input$module]]
    
    # Set the QC parameters.
    qc_params <-
      list(minSegmentReads = input$min_segment_reads, # Minimum number of reads (1000)
           percentTrimmed = input$percent_trimmed,    # Minimum % of reads trimmed (80%)
           percentStitched = input$percent_stitched,   # Minimum % of reads stitched (80%)
           percentAligned = input$percent_aligned,    # Minimum % of reads aligned (80%)
           percentSaturation = input$percent_saturation, # Minimum sequencing saturation (50%)
           minNegativeCount = input$min_negative_count,   # Minimum negative control counts (10)
           maxNTCCount = input$max_ntc_count,     # Maximum counts observed in NTC well (1000)
           minNuclei = input$min_nuclei,         # Minimum # of nuclei estimated (100)
           minArea = input$min_area)         # Minimum segment area (5000)
    data_object <-
      setSegmentQCFlags(data_object, 
                        qcCutoffs = qc_params)
    
    # Collate QC Results
    qc_results <- protocolData(data_object)[["QCFlags"]]
    flag_columns <- colnames(qc_results)
    qc_summary <- data.frame(Pass = colSums(!qc_results[, flag_columns]),
                             Warning = colSums(qc_results[, flag_columns]))
    qc_results$QCStatus <- apply(qc_results, 1L, function(x) {
      ifelse(sum(x) == 0L, "PASS", "WARNING")
    })
    qc_summary["TOTAL FLAGS", ] <-
      c(sum(qc_results[, "QCStatus"] == "PASS"),
        sum(qc_results[, "QCStatus"] == "WARNING"))
    
    
    # Get rownames.
    qc_summary <- qc_summary %>% dplyr::mutate(Parameter = rownames(.)) %>% dplyr::relocate(Parameter, .before = 1)
    
    # Convert to integers.
    qc_summary$Pass <- as.integer(qc_summary$Pass)
    qc_summary$Warning <- as.integer(qc_summary$Warning)
    
    qc_summary
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
