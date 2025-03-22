## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
target_data_object_list <- readRDS(cl_args[5])

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##
## Normalization ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Explore the relationship between the upper quartile (Q3) of the counts in each segment with the geometric mean of the negative control probes in the data. Ideally, there should be a separation between these two values to ensure we have stable measure of Q3 signal. If you do not see sufficient separation between these values, you may consider more aggressive filtering of low signal segments/genes.

### Normalization ----
# Initialize the list to hold the plots.
plot_list_normalization <- list()

for(module in names(target_data_object_list)) {
  # Check if the combined module (WTA+TCR) exists.
  # If it does, skip the individual WTA, TCR modules (because the WTA and TCR modules have to be normalized together.)
  if((combined_module %in% names(target_data_object_list)) && (module %in% c(main_module, module_tcr))) next
  
  target_data_object <- target_data_object_list[[module]]
  plot_list_normalization[[module]] <- list()
  
  # Get the negative probes. This code was used in qc_probes.R, and
  # maybe in the future, we will have a way to remove this redundant code. Maybe by saving the negative probes
  # into the NanoStringGeoMxSet object?
  negativeProbefData <- subset(fData(target_data_object), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)
  # Graph Q3 value vs negGeoMean of Negatives.
  ann_of_interest <- ann_of_interest
  if(length(neg_probes) > 1) {
    neg_probe_exprs <- t(exprs(target_data_object)[neg_probes, ])
  } else {
    neg_probe_exprs <- exprs(target_data_object)[neg_probes, ]
  }
  
  stat_data <- 
    data.frame(row.names = colnames(exprs(target_data_object)),
               Segment = colnames(exprs(target_data_object)),
               Annotation = pData(target_data_object)[, ann_of_interest],
               Q3 = unlist(apply(exprs(target_data_object), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = neg_probe_exprs) # t() because otherwise the dimensions are wrong. EDIT 2024/12/12: NO. EDIT 2025/03/14: t() when length(neg_probes) > 1, no t otherwise.
  
  stat_data_m_list <- list()
  plot_list_normalization[[module]][["Q3_norm"]] <- list()
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
      labs(x = "Counts", y = "Segments, #", title = paste0(module)) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    
    # Negative probe geometric mean (x-axis) vs Q3 value of counts of non-negative probes (y-axis). The dashed line indicates where the points (each point = 1 segment) would fall if there were no separation - i.e., if, for a given point, the Q3 values were the same as the negative probe geometric mean. Therefore, we want the points to be _above_ the dashed line.
    plt2 <- ggplot(stat_data,
                   aes(x = !!as.name(probeset), y = Q3, color = Annotation)) +
      geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
      geom_point() + guides(color = "none") + theme_bw() +
      scale_x_continuous(trans = "log2") + 
      scale_y_continuous(trans = "log2") +
      theme(aspect.ratio = 1,
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      labs(x = "Negative probe geomean, counts", y = "Q3 value, counts", title = paste0(module))
    
    plt3 <- ggplot(stat_data,
                   aes(x = !!as.name(probeset), y = Q3 / !!as.name(probeset), color = Annotation)) +
      geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
      geom_point() + theme_bw() +
      scale_x_continuous(trans = "log2") + 
      scale_y_continuous(trans = "log2") +
      theme(aspect.ratio = 1,
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      labs(x = "Negative probe geomean, counts", y = "Q3/negprobe value, counts", title = paste0(module))
    
    btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                         rel_widths = c(0.43,0.57))
    plot_list_normalization[[module]][["Q3_norm"]][[probeset]] <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
  }
  
  # Other normalization practices: https://bioconductor.org/packages/release/bioc/vignettes/GeoDiff/inst/doc/Workflow_WTA_kidney.html
  
  # Normalize.
  # https://rdrr.io/github/Nanostring-Biostats/GeomxTools/man/normalize-NanoStringGeoMxSet-method.html
  # Q3 norm (75th percentile) for WTA/CTA with or without custom spike-ins.
  target_data_object <- NanoStringNCTools::normalize(target_data_object ,
                                                     norm_method = "quant", 
                                                     desiredQuantile = .75,
                                                     toElt = "q3_norm")
  
  # Background normalization for WTA/CTA without custom spike-in.
  target_data_object <- NanoStringNCTools::normalize(target_data_object ,
                                                     norm_method = "neg",
                                                     fromElt = "exprs",
                                                     toElt = "neg_norm")
  
  # Background-subtraction correction (not used as a complete normalization method).
  target_data_object <- NanoStringNCTools::normalize(target_data_object ,
                                                     norm_method = "subtractBackground",
                                                     fromElt = "exprs",
                                                     toElt = "bg_sub")
  
  # Q3 normalization of background-subtracted data.
  target_data_object <- NanoStringNCTools::normalize(target_data_object,
                                                     norm_method = "quant",
                                                     desiredQuantile = .75,
                                                     fromElt = "bg_sub",
                                                     toElt = "bg_sub_q3")
  
  # 90th-percentile normalization of background-subtracted data (used in TCR analysis).
  target_data_object <- NanoStringNCTools::normalize(target_data_object,
                                                     norm_method = "quant",
                                                     desiredQuantile = .9,
                                                     fromElt = "bg_sub",
                                                     toElt = "bg_sub_p90")
  
  # Background normalization of background-subtracted data.
  target_data_object <- NanoStringNCTools::normalize(target_data_object,
                                                     norm_method = "neg",
                                                     fromElt = "bg_sub",
                                                     toElt = "bg_sub_neg")
  
  # Quantile normalization.
  # https://www.statology.org/quantile-normalization-in-r/
  # Quantile normalization (in which distributions are forced to be the same* across samples) is NOT the same thing as 
  # quantile-specific normalization, in which only a particular quantile (usually a quartile) is forced
  # to be the same across samples. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171491/. 
  # *tied values notwithstanding.
  quant_normalized <- normalize.quantiles(target_data_object@assayData$exprs)
  rownames(quant_normalized) <- rownames(target_data_object@assayData$exprs)
  colnames(quant_normalized) <- colnames(target_data_object@assayData$exprs)
  assayDataElement(object = target_data_object, elt = "quant", validate = FALSE) <- quant_normalized
  
  # 2024/11/15: probe removal has been moved here from qc_probes.R.
  # Subset object to exclude manually selected probes.
  probes_exclude <- probes_exclude %>% str_split(",") %>% unlist()
  if(sum(probes_exclude == "None") < length(probes_exclude)) {
    probe_qc_passed_2 <- 
      subset(target_data_object, 
             !(fData(target_data_object)[["TargetName"]] %in% probes_exclude))
    dim(probe_qc_passed_2)
    target_data_object <- probe_qc_passed_2
  }
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #
  # Save back to list.
  #
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  target_data_object_list[[module]] <- target_data_object
  
  # Visualize the first 10 segments with each normalization method, before and after normalization.
  # Set the colors.
  n_colors <- length(target_data_object@assayData)
  colors <- RColorBrewer::brewer.pal(n_colors, "Set2")
  names(colors) <- names(target_data_object@assayData)
  
  norm_methods_not_raw <- names(target_data_object@assayData) %>% .[. != "exprs"]
  for(norm_method in c("exprs", norm_methods_not_raw)) { # Ensure raw goes first. 
    plot_list_normalization[[module]][[norm_method]] <- assayDataElement(target_data_object[,1:10], elt = norm_method) %>% 
      melt %>% 
      dplyr::rename(Gene = 1, Segment = 2, Count = 3) %>%
      ggplot(aes(x = Segment, y = Count)) +
      geom_boxplot(fill = colors[names(colors)==norm_method]) +
      scale_y_continuous(trans='log10') + 
      scale_x_discrete(label = 1:10) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      labs(y = paste0("Counts, ", normalization_names[names(normalization_names)==norm_method]), x = "Segment", title = paste0(module))
  }
}

### Log-transformation ----
for(module in names(target_data_object_list)) {
  # Check if the combined module (WTA+TCR) exists.
  # If it does, skip the individual WTA, TCR modules.
  if((combined_module %in% names(target_data_object_list)) && (module %in% c(main_module, module_tcr))) next
  
  target_data_object <- target_data_object_list[[module]]
  
  for(norm_method in names(target_data_object@assayData)) {
    # Add the log-transformed data to the data object as well under the appropriate normalizations.
    assayDataElement(object = target_data_object, elt = paste0("log_", norm_method), validate = FALSE) <- # Have to set validate to FALSE; otherwise it thinks the dimensions aren't the same. ... 
      assayDataApply(target_data_object, 2, FUN = function(x) log2(x+1), elt = norm_method)
  }
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #
  # Save back to list.
  #
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  target_data_object_list[[module]] <- target_data_object
  rm(target_data_object)
  gc()
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Split combined WTA+TCR ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(combined_module %in% names(target_data_object_list)) {
  # Extract unique modules from the WTA+TCR NanoStringGeoMxSet.
  modules <- fData(target_data_object_list[[paste(c(main_module, module_tcr), collapse = ",")]]) %>% .$Module %>% unique
  
  # Split the NanoStringGeoMxSet object by Module
  individual_module_list <- lapply(modules, function(mod) {
    subset(target_data_object_list[[paste(c(main_module, module_tcr), collapse = ",")]], Module == mod)
  })
  names(individual_module_list) <- modules
  
  # Combine with `target_data_object_list`.
  target_data_object_list <- c(individual_module_list, target_data_object_list[!(names(target_data_object_list) %in% c(combined_module, modules))])
  # Removes the combined module and ensures that for the individual modules, only the normalized objects are saved.
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save the normalized NanoStringGeoMxSet to RDS.
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_normalized.rds"))
# Save the graphs.
saveRDS(plot_list_normalization, paste0(output_dir_rdata, "normalization_plot_list.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)