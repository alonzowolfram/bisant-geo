## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

# Access the PKC files, to ensure that expected PKCs have been loaded for this study.
# This code was used in qc_segments.R, and
# maybe in the future, we will have a way to remove this redundant code. Maybe by saving the modules
# into the NanoStringGeoMxSet object?
pkcs <- annotation(data_object)
modules <- base::gsub(".pkc", "", pkcs)
pkc_summary <- data.frame(PKCs = pkcs, modules = modules)

# Create the object to hold probe QC plots.
plot_list_probe_qc <- list()

# Set probe QC flags.
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers.
data_object <- setBioProbeQCFlags(data_object, 
                                  qcCutoffs = list(minProbeRatio = min_probe_ratio,
                                                   percentFailGrubbs = percent_fail_grubbs), 
                                  removeLocalOutliers = TRUE)

probe_qc_results <- fData(data_object)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(probe_qc_results[, -1]) == 0),
                    Global = sum(probe_qc_results$GlobalGrubbsOutlier),
                    Local = sum(rowSums(probe_qc_results[, -2:-1]) > 0
                                & !probe_qc_results$GlobalGrubbsOutlier))
# qc_df

# Subset object to exclude all that did not pass Ratio & Global testing.
probe_qc_passed <- 
  subset(data_object, 
         fData(data_object)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data_object)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(probe_qc_passed)
data_object <- probe_qc_passed 

# Create gene-level count data. 
# Check how many unique targets the object has.
length(unique(featureData(data_object)[["TargetName"]]))

# Collapse to targets.
target_data_object <- aggregateCounts(data_object)
dim(target_data_object)
# View the first few rows of the expression matrix.
exprs(target_data_object)[1:5, 1:2]

# Limit of quantification. 
# The LOQ is calculated based on the distribution of negative control probes and is intended to approximate the quantifiable limit of gene expression per segment. 
# Please note that this process is more stable in larger segments. 
# Likewise, the LOQ may not be as accurately reflective of true signal detection rates in segments with low negative probe counts (ex: <2).
# Define LOQ SD threshold and minimum value.
cutoff <- loq_sd_cutoff
min_loq <- min_loq

# Calculate LOQ per module tested.
loq <- data.frame(row.names = colnames(target_data_object))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_data_object)))) {
    loq[, module] <-
      pmax(min_loq,
           pData(target_data_object)[, vars[1]] * 
             pData(target_data_object)[, vars[2]] ^ cutoff)
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
# Filter out segments with exceptionally low signal ([small fraction of panel genes detected above LOQ] relative to [other segments in the study]).
# Save detection rate information to phenodata.
pData(target_data_object)$GenesDetected <- 
  colSums(loq_mat, na.rm = TRUE)
pData(target_data_object)$GeneDetectionRate <-
  pData(target_data_object)$GenesDetected / nrow(target_data_object)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%.
pData(target_data_object)$DetectionThreshold <- 
  cut(pData(target_data_object)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# Stacked bar plot of different cut points (1%, 5%, 10%, 15%).
plot_list_probe_qc[["gene_detection_rate_by_segment"]] <- ggplot(pData(target_data_object),
                                                                 aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")
#> Warning: The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
#> ℹ Please use `after_stat(count)` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.

# Create a table to review what tissue type (tumor vs normal) is going to be impacted by each threshold:
# Cut percent genes detected at 1, 5, 10, 15.
kable(table(pData(target_data_object)$DetectionThreshold,
            pData(target_data_object)$Segment))

# Remove segments with <5% of genes detected.
target_data_object <-
  target_data_object[, pData(target_data_object)$GeneDetectionRate >= .05]

dim(target_data_object)

# Re-plot the Sankey diagram showing our current working dataset. This is now a dataset that no longer contains segments flagged by Segment QC or that have low gene detection rates.
# 
# # select the annotations we want to show, use `` to surround column names with
# # spaces or special symbols
# count_mat <- count(pData(data_object), `slide name`, class, region, segment)
# # simplify the slide names
# count_mat$`slide name` <- 
#     gsub("disease", "d",
#          gsub("normal", "n", count_mat$`slide name`))
# # gather the data and plot in order: class, slide name, region, segment
# test_gr <- gather_set_data(count_mat, 1:4)
# test_gr$x <-
#     factor(test_gr$x,
#            levels = c("class", "slide name", "region", "segment"))
# # plot Sankey
# ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
#     geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
#     geom_parallel_sets_axes(axis.width = 0.2) +
#     geom_parallel_sets_labels(color = "white", size = 5) +
#     theme_classic(base_size = 17) + 
#     theme(legend.position = "bottom",
#           axis.ticks.y = element_blank(),
#           axis.line = element_blank(),
#           axis.text.y = element_blank()) +
#     scale_y_continuous(expand = expansion(0)) + 
#     scale_x_discrete(expand = expansion(0)) +
#     labs(x = "", y = "") +
#     annotate(geom = "segment", x = 4.25, xend = 4.25, y = 20, 
#              yend = 120, lwd = 2) +
#     annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
#              hjust = 0.5, label = "100 segments")

# Gene detection rate. 
# Calculate detection rate:
loq_mat <- loq_mat[, colnames(target_data_object)]
fData(target_data_object)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
fData(target_data_object)$DetectionRate <-
  fData(target_data_object)$DetectedSegments / nrow(pData(target_data_object))

# # Gene of interest detection table
# goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
#          "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
# goi_df <- data.frame(
#     Gene = goi,
#     Number = fData(target_data_object)[goi, "DetectedSegments"],
#     DetectionRate = percent(fData(target_data_object)[goi, "DetectionRate"]))

# Gene filtering. 
# Graph the total number of genes detected in different percentages of segments. Based on the visualization below, we can better understand global gene detection in our study and select how many low detected genes to filter out of the dataset. Gene filtering increases performance of downstream statistical tests and improves interpretation of true biological signal.
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_data_object)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_data_object))
rownames(plot_detect) <- plot_detect$Freq

plot_list_probe_qc[["percent_segments_by_gene_detection_cutoff"]] <-ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
# Also manually include the negative control probe, for downstream use.
negativeProbefData <- subset(fData(target_data_object), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_data_object <- 
  target_data_object[fData(target_data_object)$DetectionRate >= 0.1 |
                       fData(target_data_object)$TargetName %in% neg_probes, ]
dim(target_data_object)

# # Retain only detected genes of interest.
# goi <- goi[goi %in% rownames(target_data_object)]

# Save the QC-ed NanoStringGeoMxSet object.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-probes.rds"))

# Add everything to the PowerPoint. 
# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Probe QC"), 
                   location = ph_location_label(ph_label = "Title 1"))
# Add the graphs.
for(item in names(plot_list_probe_qc)) {
  if(item != "neg_geo_means") {
    # Add to the PowerPoint. 
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Probe QC"),
                       location = ph_location_label(ph_label = "Title 1")) %>% 
      officer::ph_with(value = plot_list_probe_qc[[item]],
                       location = ph_location_label(ph_label = "Content Placeholder 2"))
  } else {
    for(ann in names(plot_list_probe_qc[["neg_geo_means"]])) {
      # Add to the PowerPoint. 
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Probe QC"),
                         location = ph_location_label(ph_label = "Title 1")) %>% 
        officer::ph_with(value = plot_list_probe_qc[[item]][[ann]],
                         location = ph_location_label(ph_label = "Content Placeholder 2"))
    }
  }
}
# # Add the summary table. 
# qc_summary <- qc_summary %>% 
#   dplyr::mutate(Metric = rownames(.)) %>% 
#   dplyr::relocate(Metric, .before = 1)
# pptx <- pptx %>% 
#   officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
#   officer::ph_with(value = paste0("Probe QC summary"),
#                    location = ph_location_label(ph_label = "Title 1")) %>% 
#   officer::ph_with(value = qc_summary,
#                    location = ph_location_label(ph_label = "Content Placeholder 2"))
# print(pptx, ppt_output_path)
print(pptx, cl_args[5])