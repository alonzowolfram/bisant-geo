## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

# Explore the relationship between the upper quartile (Q3) of the counts in each segment with the geometric mean of the negative control probes in the data. Ideally, there should be a separation between these two values to ensure we have stable measure of Q3 signal. If you do not see sufficient separation between these values, you may consider more aggressive filtering of low signal segments/genes.

# Initialize the list to hold the plots.
plot_list_normalization <- list()

# Get the negative probes. This code was used in qc_probes.R, and
# maybe in the future, we will have a way to remove this redundant code. Maybe by saving the negative probes
# into the NanoStringGeoMxSet object?
negativeProbefData <- subset(fData(target_data_object), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
# Graph Q3 value vs negGeoMean of Negatives.
ann_of_interest <- ann_of_interest
stat_data <- 
  data.frame(row.names = colnames(exprs(target_data_object)),
             Segment = colnames(exprs(target_data_object)),
             Annotation = pData(target_data_object)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_data_object), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = t(exprs(target_data_object)[neg_probes, ])) # t() because otherwise the dimensions are wrong. 

stat_data_m_list <- list()
plot_list_normalization[["Q3_norm"]] <- list()
# One entry for each negative probe set.
for(probeset in colnames(stat_data) %>% regexPipes::grep("NegProbe", value=T)) {
  stat_data_m <- melt(stat_data, measure.vars = c("Q3", probeset),
                      variable.name = "Statistic", value.name = "Value")
  stat_data_m_list[[probeset]] <- stat_data_m
  
  plt1 <- ggplot(stat_data_m,
                 aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")
  
  # Negative probe geometric mean (x-axis) vs Q3 value of counts of non-negative probes (y-axis). The dashed line indicates where the points (each point = 1 segment) would fall if there were no separation - i.e., if, for a given point, the Q3 values were the same as the negative probe geometric mean. Therefore, we want the points to be _above_ the dashed line.
  plt2 <- ggplot(stat_data,
                 aes(x = !!as.name(probeset), y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
  
  plt3 <- ggplot(stat_data,
                 aes(x = !!as.name(probeset), y = Q3 / !!as.name(probeset), color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
  
  btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                       rel_widths = c(0.43,0.57))
  plot_list_normalization[["Q3_norm"]][[probeset]] <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
}

# Normalize.
# Q3 norm (75th percentile) for WTA/CTA with or without custom spike-ins.
target_data_object <- NanoStringNCTools::normalize(target_data_object ,
                                                   norm_method = "quant", 
                                                   desiredQuantile = .75,
                                                   toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in.
target_data_object <- NanoStringNCTools::normalize(target_data_object ,
                                                   norm_method = "neg",
                                                   fromElt = "exprs",
                                                   toElt = "neg_norm")

# Visualize the first 10 segments with each normalization method, before and after normalization.
# Raw.
plot_list_normalization[["before_norm"]] <- exprs(target_data_object)[,1:10] %>% 
  melt %>% 
  dplyr::rename(Gene = 1, Segment = 2, Count = 3) %>%
  ggplot(aes(x = Segment, y = Count)) +
  geom_boxplot(fill = "#9EDAE5") +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(label = 1:10) +
  theme_bw() +
  labs(y = "Counts, Raw", x = "Segment")

# Q3 norm counts.
plot_list_normalization[["after_Q3_norm"]] <- assayDataElement(target_data_object[,1:10], elt = "q_norm") %>% 
  melt %>% 
  dplyr::rename(Gene = 1, Segment = 2, Count = 3) %>%
  ggplot(aes(x = Segment, y = Count)) +
  geom_boxplot(fill = "#2CA02C") +
  scale_y_continuous(trans='log10') + 
  scale_x_discrete(label = 1:10) +
  theme_bw() +
  labs(y = "Counts, Q3 normalized", x = "Segment")
# Negative norm counts.
plot_list_normalization[["after_neg_norm"]] <- assayDataElement(target_data_object[,1:10], elt = "neg_norm") %>% 
  melt %>% 
  dplyr::rename(Gene = 1, Segment = 2, Count = 3) %>%
  ggplot(aes(x = Segment, y = Count)) +
  geom_boxplot(fill = "#FF7F0E") +
  scale_y_continuous(trans='log10') + 
  scale_x_discrete(label = 1:10) +
  theme_bw() +
  labs(y = "Counts, negative normalized", x = "Segment")

# Add everything to the PowerPoint. 
# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Normalization"), 
                   location = ph_location_label(ph_label = "Title 1"))

# Graphing parameters.
plot_width_q3 <- 12
plot_height_q3 <- 12
plot_width <- 6
plot_height <- 6
units <- "in"
res <- 300
# Add the graphs.
for(item in names(plot_list_normalization)) {
  if(item=="Q3_norm") {
    for(item2 in names(plot_list_normalization[[item]])) {
      plot <- plot_list_normalization[[item]][[item2]]
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, "normalization_", item, "-", item2, ".eps")
      png_path <- paste0(output_dir_pubs, "normalization_", item, "-", item2, ".png")
      saveEPS(plot, eps_path, width = plot_width_q3, height = plot_height_q3)
      savePNG(plot, png_path, width = plot_width_q3, height = plot_height_q3, units = units, res = res)
      
      # Add to the PowerPoint. 
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Normalization"),
                         location = ph_location_label(ph_label = "Title 1")) %>% 
        officer::ph_with(value = external_img(png_path, width = plot_width_q3, height = plot_height_q3, unit = units),
                         location = ph_location_label(ph_label = "Content Placeholder 2"),
                         use_loc_size = FALSE)
      
    } 
  } else {
    plot <- plot_list_normalization[[item]]
    
    # Save to EPS and PNG and then ...
    eps_path <- paste0(output_dir_pubs, "normalization_", item, ".eps")
    png_path <- paste0(output_dir_pubs, "normalization_", item, ".png")
    saveEPS(plot, eps_path, width = plot_width, height = plot_height)
    savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
    
    # Add to the PowerPoint. 
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Normalization"),
                       location = ph_location_label(ph_label = "Title 1")) %>% 
      officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                       location = ph_location_label(ph_label = "Content Placeholder 2"),
                       use_loc_size = FALSE)
  }
}

## ---------------------------
# Export to disk.

# Save the normalized NanoStringGeoMxSet to RDS.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_normalized.rds"))
# Output everything to the PowerPoint. 
print(pptx, cl_args[5])