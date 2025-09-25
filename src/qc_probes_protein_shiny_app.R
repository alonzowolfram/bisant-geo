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

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Load in the data object list.
data_object_list <- readRDS("Rdata/NanoStringGeoMxSet_qc-segments.rds")
# Get the data object
data_object <- data_object_list[[1]]
exprs <- data_object@assayData$exprs

# Pull IgGs (negative controls) from data set
igg_names <- iggNames(data_object)

# Get geomean of negative controls
gmean_nc <- (exprs[igg_names,] + 1) %>% log() %>% mean(na.rm = T) %>% exp # https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

# Divide all counts by geomean of negative controls
snr <- (exprs + 1) / gmean_nc
# Then take the log
log2_snr <- log2(snr)

# Order proteins by average log2(SNR), ascending
log2_snr <- rbind(log2_snr%>% .[rownames(.) %in% igg_names,,drop=F],
                  log2_snr %>% .[!(rownames(.) %in% igg_names),,drop=F] %>% .[order(rowMeans(.), decreasing = F),,drop=F]
)
protnames <- rownames(log2_snr)

# Reshape log2_snr to long format
long <- log2_snr %>% 
  as.data.frame %>% # Convert to data frame so we can reshape it long
  tibble::rownames_to_column("Protein") %>% 
  tidyr::pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Value") %>%
  mutate(Protein = factor(Protein, levels = protnames))

# Set the number of proteins per page in graph
PAGE_SIZE <- 64

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## UI ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("QC: Probes"),
  
  # Number of proteins per page
  fluidRow(
    column(6,
           sliderInput("pagesize", "Proteins per page", min=16, max=128, value=PAGE_SIZE, step=16)
    )
  ),
  # Plot
  fluidRow(plotOutput("signal_dist")),
  # Page navigation
  fluidRow(
    column(6,
           actionButton("prev_page", "◀ Prev"),
           actionButton("next_page", "Next ▶"),
           span(style="margin-left:10px;", textOutput("page_lab", inline = TRUE))
    )
  )
  # ,
  # 
  # # Row 2: input for genes of interest
  # # Text input for genes to be highlighted
  # textAreaInput("highlight_genes",
  #               "Genes of interest to highlight on the plot"),
  # actionButton("submit_highlight_genes", "Submit")
)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Server ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# https://bioconductor.org/packages/devel/bioc/vignettes/GeomxTools/inst/doc/Protein_in_GeomxTools.html
server <- function(input, output) {
  # Pagination logic:
  # Calculate the total number of pages
  # 
  
  # Total pages depends on page size
  total_pages <- reactive({
    n <- length(levels(long$Protein))
    psize <- input$pagesize
    max(1, ceiling(n / psize))
  })
  
  # Current page index
  page <- reactiveVal(1)
  
  observeEvent(input$next_page, page(min(page() + 1, total_pages())))
  observeEvent(input$prev_page, page(max(page() - 1, 1)))
  observeEvent(input$pagesize, {
    # keep the first gene visible after page-size change
    page(1)
  })
  
  # Print current page
  output$page_lab <- renderText(sprintf("Page %d of %d", page(), total_pages()))
  
  # Slice the factor levels for the current page
  current_levels <- reactive({
    lv <- levels(long$Protein)
    start <- (page() - 1) * input$pagesize + 1
    end   <- min(start + input$pagesize - 1, length(lv))
    lv[start:end]
  })
  
  # Plot
  output$signal_dist <- renderPlot({
    # Create the plot
    # Based on qcProteinSignal(object = data_object, neg.names = igg_names)
    
    # Get proteins of interest
    protein_order <- qcProteinSignalNames(object = data_object, neg.names = igg_names)
    genes_of_interest_input <- input$highlight_genes %>% str_split(",") %>% unlist %>% intersect(protein_order)
    genes_of_interest <- c()
    if(length(genes_of_interest_input) > 0) {
      for(gene in genes_of_interest_input) {
        genes_of_interest <- c(genes_of_interest, which(protein_order==gene))
      } 
    }
    
    # Display the plot with genes of interest highlighted
    long %>%
      dplyr::filter(Protein %in% current_levels()) %>%
      dplyr::mutate(Protein = factor(Protein, levels = current_levels())) %>%
      ggplot(aes(x = Protein, y = Value)) +
      geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +  # boxplots
      geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "blue") +    # jittered points
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      xlab("") + 
      ylab("Log2 signal-to-background ratio")
    # for(i in 1:length(genes_of_interest)) {
    #   xleft <- ifelse(i==1, 0, genes_of_interest[i]-1)
    #   rect()
    # }
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
