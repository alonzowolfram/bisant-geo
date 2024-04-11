## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
data_object <- readRDS(cl_args[4])
pkcs <- annotation(data_object)
modules <- base::gsub(".pkc", "", pkcs)
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

## ---------------------------

# QC: segments
# Select the QC parameter cutoffs, against which ROI/AOI segments will be tested and flagged appropriately.
qc_params <-
  list(minSegmentReads = min_segment_reads, # Minimum number of reads (1000)
       percentTrimmed = percent_trimmed,    # Minimum % of reads trimmed (80%)
       percentStitched = percent_stitched,   # Minimum % of reads stitched (80%)
       percentAligned = percent_aligned,    # Minimum % of reads aligned (80%)
       percentSaturation = percent_saturation, # Minimum sequencing saturation (50%)
       minNegativeCount = min_negative_count,   # Minimum negative control counts (10)
       maxNTCCount = max_ntc_count,     # Maximum counts observed in NTC well (1000)
       minNuclei = min_nuclei,         # Minimum # of nuclei estimated (100)
       minArea = min_area)         # Minimum segment area (5000)
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

# Visualize segment QC.
col_by <- "segment"

# Graphical summaries of QC statistics plot function
# Generate the QC histograms.
plot_list_segment_qc <- list()
# Trimmed %
plot_list_segment_qc[["trimmed"]] <- QCHistogram(sData(data_object), "Trimmed (%)", col_by, percent_trimmed)
# Stitched %
plot_list_segment_qc[["stitched"]] <- QCHistogram(sData(data_object), "Stitched (%)", col_by, percent_stitched)
# Aligned %
plot_list_segment_qc[["aligned"]] <- QCHistogram(sData(data_object), "Aligned (%)", col_by, percent_aligned)
# Sequencing saturation %
plot_list_segment_qc[["seq_sat"]] <- QCHistogram(sData(data_object), "Saturated (%)", col_by, percent_saturation) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
# Area.
plot_list_segment_qc[["area"]] <- QCHistogram(sData(data_object), "area", col_by, min_area, scale_trans = "log10")
# Nuclei.
plot_list_segment_qc[["nuclei"]] <- QCHistogram(sData(data_object), "nuclei", col_by, min_nuclei)
# Negative geometric means.
plot_list_segment_qc[["neg_geo_means"]] <- list()
# Calculate the negative geometric means for each module.
negativeGeoMeans <- 
  esBy(negativeControlSubset(data_object), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(data_object)[["NegGeoMean"]] <- negativeGeoMeans
# Explicitly copy the Negative geoMeans from sData to pData.
negCols <- paste0("NegGeoMean_", modules)
pData(data_object)[, negCols] <- sData(data_object)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QCHistogram(pData(data_object), ann, col_by, min_negative_count, scale_trans = "log10")
  plot_list_segment_qc[["neg_geo_means"]][[ann]] <- plt 
  # print(plt)
}
# If there are NTCs in the experiment (sData(data_object)$NTC does not return NULL), add NTC stats.
# # No template control (NTC)
# # Detatch neg_geomean columns ahead of aggregateCounts call.
if(!is.null(sData(data_object)$NTC)) {
  pData(data_object) <- pData(data_object)[, !colnames(pData(data_object)) %in% negCols]

  # Show all NTC values, Freq = # of Segments with a given NTC count:
  ntc_table <- as.data.frame(table(NTC_Count = sData(data_object)$NTC),
      col.names = c("NTC Count", "# of Segments"))
}
# Summarize all QC information in a table.
qc_summary_kable <- kable(qc_summary, caption = "QC Summary Table for each Segment")

# Remove flagged segments that do not meet QC cutoffs.
data_object <- data_object[, qc_results$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(data_object)

# Add everything to the PowerPoint. 
# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Segment QC"), 
                   location = ph_location_label(ph_label = "Title 1"))

# Graphing parameters.
plot_width <- 7
plot_height <- 7 
units <- "in"
res <- 300
# Add the graphs.
for(item in names(plot_list_segment_qc)) {
  if(item != "neg_geo_means") {
    plot <- plot_list_segment_qc[[item]]
    
    # Save to EPS and PNG and then ...
    eps_path <- paste0(output_dir_pubs, "qc-segments_", item, ".eps")
    png_path <- paste0(output_dir_pubs, "qc-segments_", item, ".png")
    saveEPS(plot, eps_path, width = plot_width, height = plot_height)
    savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
    
    # Add to the PowerPoint. 
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Segment QC"),
                       location = ph_location_label(ph_label = "Title 1")) %>% 
      officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                       location = ph_location_label(ph_label = "Content Placeholder 2"),
                       use_loc_size = FALSE) # use_loc_size = FALSE forces the size to be the size of the saved PNG, not the PPT slide. 
  } else {
    for(ann in names(plot_list_segment_qc[["neg_geo_means"]])) {
      plot <- plot_list_segment_qc[[item]][[ann]]
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, "qc-segments_", item, "-", ann, ".eps")
      png_path <- paste0(output_dir_pubs, "qc-segments_", item, "-", ann, ".png")
      saveEPS(plot, eps_path, width = plot_width, height = plot_height)
      savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
      
      # Add to the PowerPoint. 
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Segment QC"),
                         location = ph_location_label(ph_label = "Title 1")) %>% 
        officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                         location = ph_location_label(ph_label = "Content Placeholder 2"),
                         use_loc_size = FALSE) # use_loc_size = FALSE forces the size to be the size of the saved PNG, not the PPT slide. 
    }
  }
}
# Add the NTC table if it exists.
if(exists("ntc_table")) {
  pptx <- pptx %>% 
    officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
    officer::ph_with(value = paste0("NTC summary"),
                     location = ph_location_label(ph_label = "Title 1")) %>% 
    officer::ph_with(value = ntc_table,
                     location = ph_location_label(ph_label = "Content Placeholder 2"))
}
# Add the summary table. 
qc_summary <- qc_summary %>% 
  dplyr::mutate(Metric = rownames(.)) %>% 
  dplyr::relocate(Metric, .before = 1)
pptx <- pptx %>% 
  officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Segment QC summary"),
                   location = ph_location_label(ph_label = "Title 1")) %>% 
  officer::ph_with(value = qc_summary,
                   location = ph_location_label(ph_label = "Content Placeholder 2"))

## ---------------------------
# Export to disk.

# Save the QC-ed NanoStringGeoMxSet object.
saveRDS(data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-segments.rds"))
# Output everything to the PowerPoint. 
print(pptx, cl_args[5])