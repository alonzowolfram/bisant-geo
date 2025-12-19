## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
data_object_list <- readRDS(cl_args[5])

# Access the PKC files, to ensure that expected PKCs have been loaded for this study
# This code was used in qc_segments.R, and
# maybe in the future, we will have a way to remove this redundant code. Maybe by saving the modules
# into the NanoStringGeoMxSet object?
modules <- names(data_object_list)
pkcs <- paste0(modules, ".pkc")
pkc_summary <- data.frame(PKCs = pkcs, modules = modules)

# Set `main_module` if not set already
if(flagVariable(main_module)) main_module <- modules[1]

# Set the probes to be included no matter what
probes_include <- probes_include %>% str_split(",") %>% unlist()

# Set up `target_data_object_list`, `data_object_fData`, `plot_list_probe_qc`, and `qc_df` for protein
# We will only overwrite these (and perform probe QC) if the analyte is RNA
target_data_object_list <- data_object_list
data_object_fData <- fData(data_object_list[[1]])
plot_list_probe_qc <- list()
qc_df <- data.frame()

if(analyte=="protein") {
  target_data_object <- target_data_object_list[[main_module]]
    
  # Get the IgG names, which we will use to determine the background expression level
  igg_names <- iggNames(target_data_object)
  hk_names <- hkNames(target_data_object)
  # Extract the protein expression matrix
  exp_protein <- target_data_object@assayData$exprs
  
  # Perform probe QC only if `percent_of_segments` is set and numeric
  if(is.numeric(percent_of_segments)) {
    # Get the per-sample background
    sample_bg <- exp_protein %>% 
      t %>% 
      as.data.frame() %>% 
      mutate(IgG_bg = exp(rowMeans(log(select(., all_of(igg_names)) + 1)))) %>% 
      tibble::rownames_to_column(var = "sample") %>% 
      dplyr::mutate(sample = sample %>% regexPipes::gsub("-[[:alpha:]]{1}-", "-")) %>% 
      dplyr::select(sample, IgG_bg)
    
    # Calculate signal:noise ratios
    df <- exp_protein %>% t %>% as.data.frame %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "sample") %>% 
      dplyr::mutate(sample = sample %>% regexPipes::gsub("-[[:alpha:]]{1}-", "-")) %>%
      # Get into long form as required by ggplot2
      reshape2::melt(id.vars = c("sample"), variable.name = "protein", value.name = "expression") %>%
      # Left join BG
      dplyr::left_join(sample_bg) %>% 
      # Get the log2 of the signal:BG ratio
      dplyr::mutate(log2SignalBgRatio = log2(expression + 1) - log2(IgG_bg + 1))
    # Make a list of proteins with signal:noise ratio >= 2 in < `percent_of_segments`% of ROIs
    n_roi_total <- nrow(pData(target_data_object))
    protein_qc_table <- df %>% 
      group_by(protein) %>% 
      dplyr::filter(log2SignalBgRatio >= min_log2_protein_snr) %>% 
      dplyr::summarise(n_pass = n()) %>% 
      dplyr::mutate(percent_pass = n_pass / !!n_roi_total * 100, MinLog2SignalNoiseRatio = min_log2_protein_snr)
    proteins_fail <- protein_qc_table %>% dplyr::filter(percent_pass < !!percent_of_segments) %>% .[["protein"]] %>% as.character
    proteins_fail <- setdiff(proteins_fail, union(igg_names, hk_names)) %>% setdiff(probes_include) # Also keep any probes specified by `probes_include`
    # View(proteins_fail %>% as.data.frame)
    
    # Remove proteins that fail and are not included in `probes_include`
    # https://bioconductor.org/packages/release/bioc/vignettes/GeomxTools/inst/doc/Developer_Introduction_to_the_NanoStringGeoMxSet.html
    # dim(subset(target_data_object, !(TargetName %in% proteins_fail)))
    target_data_object <- subset(target_data_object, !(TargetName %in% proteins_fail))
    
    # Save the QC-ed target data object back to the list
    target_data_object_list[[main_module]] <- target_data_object
    
    # Create a histogram of protein signal:noise ratio
    # Custom breakpoints
    breaks <- c(0, 5, 10, 25, 50, 75, 100)
    protein_qc_table$bin <- cut(protein_qc_table$percent_pass, breaks = breaks, right = F)
    protein_qc_table <- protein_qc_table %>% 
      dplyr::mutate(bin = factor(
        case_match(
          bin,
          "[0,5)" ~ "0 - 5%",
          "[5,10)" ~ "5 - 10%",
          "[10,25)" ~ "10 - 25%",
          "[25,50)" ~ "25 - 50%",
          "[50,75)" ~ "50 - 75%",
          "[75,100)" ~ "75 - 100%"),
        levels = c("0 - 5%", "5 - 10%", "10 - 25%", "25 - 50%", "50 - 75%", "75 - 100%")
      )) %>%
      dplyr::mutate(bin = replace_na(bin, "75 - 100%"),
                    Class = "All proteins")
      
    # Graph
    plot_list_probe_qc[["protein_detection_rate_by_segment"]] <- protein_qc_table %>% ggplot(aes(x = bin)) +
      geom_bar(aes(fill = Class)) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(y = "# of proteins", x = glue::glue("% ROIs passed (log2(signal:noise ratio) >= {min_log2_protein_snr})"))
    
  }
  
} else if(analyte=="RNA") {
  ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ##                                                                
  ## Probe QC ----
  ##
  ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Create a list to hold the target data objects.
  target_data_object_list <- list()
  # Create a list to hold the probe QC table.
  qc_df_list <- list()
  
  # Loop over each module and calculate probe QC metrics.
  if(flagVariable(modules_filter_probes)) modules_filter_probes <- c("")
  
  for(module in names(data_object_list)) {
    data_object <- data_object_list[[module]]
    
    # Set probe QC flags.
    # Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
    # FALSE if you do not want to remove local outliers.
    data_object <- setBioProbeQCFlags(data_object, 
                                      qcCutoffs = list(minProbeRatio = min_probe_ratio,
                                                       percentFailGrubbs = percent_fail_grubbs), 
                                      removeLocalOutliers = TRUE)
    
    probe_qc_results <- fData(data_object)[["QCFlags"]]
    data_object_fData <- fData(data_object)
    
    # Define QC table for Probe QC
    qc_df <- data.frame(Module = module,
                        Passed = sum(rowSums(probe_qc_results[, -1]) == 0),
                        Global = sum(probe_qc_results$GlobalGrubbsOutlier),
                        Local = sum(rowSums(probe_qc_results[, -2:-1]) > 0
                                    & !probe_qc_results$GlobalGrubbsOutlier))
    qc_df_list[[module]] <- qc_df
    
    # If module is in modules_filter_probes, 
    # subset object to exclude all that did not pass Ratio & Global testing
    # But keep all the ones we specifically want to include, including TCR and negative probes
    if(module %in% modules_filter_probes) {
      negativeProbefData <- subset(fData(data_object), CodeClass == "Negative")
      neg_probes <- unique(negativeProbefData$TargetName)
      probe_qc_passed <- 
        subset(data_object, 
               (fData(data_object)[["QCFlags"]][,c("LowProbeRatio")] == FALSE & fData(data_object)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE) | # Include probes that passed QC.
                 fData(data_object)$TargetName %in% c(probes_include, neg_probes) | # Include specified probes and negative probes
                 base::grepl("TR[A/B/D/G][C/J/V]", fData(data_object)$TargetName) # Include TCR probes
        )
      dim(probe_qc_passed)
      data_object <- probe_qc_passed
    }
    
    # 2024/11/15: manual probe exclusion has been moved to normalization.R
    # # Subset object to exclude manually selected probes
    # probes_exclude <- probes_exclude %>% str_split(",") %>% unlist()
    # if(sum(probes_exclude == "None") < length(probes_exclude)) {
    #   probe_qc_passed_2 <- 
    #     subset(data_object, 
    #            !(fData(data_object)[["TargetName"]] %in% probes_exclude))
    #   dim(probe_qc_passed_2)
    #   data_object <- probe_qc_passed_2
    # }
    
    # Create gene-level count data
    # Check how many unique targets the object has
    length(unique(featureData(data_object)[["TargetName"]]))
    
    # Collapse to targets
    target_data_object <- aggregateCounts(data_object)
    dim(target_data_object)
    # View the first few rows of the expression matrix
    exprs(target_data_object)[1:5, 1:2]
    
    # Save the updated data object back to the list
    target_data_object_list[[module]] <- target_data_object
    
    # Clean up
    rm(data_object, target_data_object)
    gc()
  }
  # Merge the QC tables into one
  qc_df <- do.call(rbind, qc_df_list)
  
  ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ##                                                                
  ## Target QC ----
  ##
  ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # Limit of quantification
  # The LOQ is calculated based on the distribution of negative control probes and is intended to approximate the quantifiable limit of gene expression per segment
  # Please note that this process is more stable in larger segments
  # Likewise, the LOQ may not be as accurately reflective of true signal detection rates in segments with low negative probe counts (ex: <2)
  # Define LOQ SD threshold and minimum value
  cutoff <- loq_sd_cutoff
  min_loq <- min_loq
  
  # Calculate LOQ per module tested
  # However, for the combined module(s), if it exists, we will just use the WTA (main) module's LOQ
  for(module in modules) {
    if(base::grepl("combined_module", module)) next
    
    target_data_object <- target_data_object_list[[module]]
    
    loq <- data.frame(row.names = colnames(target_data_object))
    
    vars <- expand.grid(V1 = c("NegGeoMean", "NegGeoSD"),
                        V2 = (module %>% str_split(",") %>% unlist)) %>% dplyr::mutate(vars = paste(V1, V2, sep = "_")) %>% .$vars
    if(all(vars[1:2] %in% colnames(pData(target_data_object)))) {
      loq[, module] <-
        pmax(min_loq,
             pData(target_data_object)[, vars[1]] * 
               pData(target_data_object)[, vars[2]] ^ cutoff)
    }
    
    pData(target_data_object)$LOQ <- loq
    
    # Add back to the list
    target_data_object_list[[module]] <- target_data_object
    
    # Clean up
    rm(target_data_object)
    gc()
  }
  # If any combined modules exist, use the LOQ from the WTA module
  combined_modules <- names(target_data_object_list) %>% regexPipes::grep("combined_module", value = TRUE)
  if(length(combined_modules) > 0) {
    for(combined_module in combined_modules) {
      target_data_object <- target_data_object_list[[combined_module]]
      
      pData(target_data_object)$LOQ <- pData(target_data_object_list[[main_module]])$LOQ
      
      # Add back to the list
      target_data_object_list[[combined_module]] <- target_data_object
      
      # Clean up
      rm(target_data_object)
      gc()
    }
  }
  
  # Calculate detection rates
  loq_mat_list <- list()
  goi_df_list <- list()
  for(module in modules) {
    # Again, for the combined module, we will just use the WTA QC metrics
    if(base::grepl("combined_module", module)) next
    
    target_data_object <- target_data_object_list[[module]]
    loq <- pData(target_data_object)$LOQ
    
    loq_mat <- c()
    ind <- fData(target_data_object)$Module %in% (module %>% str_split(",") %>% unlist) # fData(target_data_object)$Module == module
    Mat_i <- t(esApply(target_data_object[ind, ], MARGIN = 1,
                       FUN = function(x) {
                         x > loq[, module]
                       }))
    loq_mat <- rbind(loq_mat, Mat_i)
    # Ensure ordering since this is stored outside of the geomxSet
    loq_mat <- loq_mat[fData(target_data_object)$TargetName, ]
    
    # Save to the list
    loq_mat_list[[module]] <- loq_mat
    
    # Save detection rate information to phenodata
    pData(target_data_object)$GenesDetected <- 
      colSums(loq_mat, na.rm = TRUE)
    pData(target_data_object)$GeneDetectionRate <-
      pData(target_data_object)$GenesDetected / nrow(target_data_object)
    
    # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
    pData(target_data_object)$DetectionThreshold <- 
      cut(pData(target_data_object)$GeneDetectionRate,
          breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
          labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
    
    # Gene detection rate
    # Calculate detection rate:
    fData(target_data_object)$DetectedSegments <- rowSums(loq_mat, na.rm = TRUE)
    fData(target_data_object)$DetectionRate <-
      fData(target_data_object)$DetectedSegments / nrow(pData(target_data_object))
    
    # Gene of interest detection table
    if(!flagVariable(genes_of_interest)) {
      goi_df <- data.frame(
        Gene = fData(target_data_object)[genes_of_interest, "TargetName"],
        Number = fData(target_data_object)[genes_of_interest, "DetectedSegments"],
        DetectionRate = percent(fData(target_data_object)[genes_of_interest, "DetectionRate"]))
      goi_df_list[[module]] <- goi_df
      # If any genes in `genes_of_interest` are not present in `rownames(fData(target_data_object))`, they will appear as `NA`s
    }
    
    # Save the target_data_object back to the list
    target_data_object_list[[module]] <- target_data_object
    
    # Clean up
    rm(target_data_object)
    if(exists("goi_df")) rm(goi_df)
    gc()
  }
  if(!flagVariable(genes_of_interest)) {
    # Check if `goi_df_list` has any (non-NA) elements
    goi_df <- do.call(cbind, goi_df_list)
    if(!is.null(goi_df)) { goi_df <- goi_df %>% dplyr::filter(if_any(everything(), ~ !is.na(.)))  } # https://stackoverflow.com/a/50131415 
    if(is.null(goi_df) || nrow(goi_df) < 1) {
      warning("None of the genes of interest you specified in `genes_of_interest` were found in the featureData. No genes of interest data frame will be returned")
      rm(goi_df)
    } 
  }
  # If any combined modules exist, use the detection rates from the WTA module
  if(length(combined_modules) > 0) {
    for(combined_module in combined_modules) {
      # Get the columns that are in the pData for the main module but missing in the pData for the combined module
      cols <- setdiff(colnames(pData(target_data_object_list[[main_module]])), colnames(pData(target_data_object_list[[combined_module]])))
      # Add them to the combined module
      pData(target_data_object_list[[combined_module]]) <- cbind(pData(target_data_object_list[[combined_module]]), pData(target_data_object_list[[main_module]])[,cols])
      fData(target_data_object_list[[combined_module]]) <- fData(target_data_object_list[[combined_module]]) %>% 
        dplyr::left_join(fData(target_data_object_list[[main_module]]) %>% 
                           dplyr::select(TargetName, DetectedSegments, DetectionRate), 
                         by = "TargetName")
    }
  }
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #
  # For the remaining part, target QC will be performed only on the main (transcriptomic) module
  #
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Get the target_data_object for the main module
  # if(combined_module %in% names(target_data_object_list)) {
  #   module <- combined_module
  # } else {
  #   module <- main_module
  # }
  module <- main_module
  target_data_object <- target_data_object_list[[module]]
  
  # Create the object to hold probe QC plots
  plot_list_probe_qc <- list()
  
  ## ................................................
  ##
  ### Segment filtering ----
  ##
  ## ................................................
  # Filter out segments with exceptionally low signal ([small fraction of panel genes detected above LOQ] relative to [other segments in the study])
  
  # Stacked bar plot of different cut points (1%, 5%, 10%, 15%)
  plot_list_probe_qc[["gene_detection_rate_by_segment"]] <- ggplot(pData(target_data_object),
                                                                   aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = segment)) +
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
  # Cut percent genes detected at 1, 5, 10, 15
  kable(table(pData(target_data_object)$DetectionThreshold,
              pData(target_data_object)$segment))
  
  # Remove segments with <x% of genes detected
  segments_keep <- pData(target_data_object)$GeneDetectionRate >= (gene_detection_rate/100) # We will later use this to subset the target data objects for the other modules
  target_data_object <-
    target_data_object[, segments_keep] # gene_detection_rate is given as a percentage, not a decimal, so we need to divide by 100 first
  dim(target_data_object)
  
  ## ................................................
  ##
  ### Gene filtering ----
  ##
  ## ................................................
  # Graph the total number of genes detected in different percentages of segments. Based on the visualization below, we can better understand global gene detection in our study and select how many low detected genes to filter out of the dataset. Gene filtering increases performance of downstream statistical tests and improves interpretation of true biological signal
  # Plot detection rate:
  plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
  plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_data_object)$DetectionRate >= x, na.rm = TRUE)}))
  plot_detect$Rate <- plot_detect$Number / nrow(fData(target_data_object))
  rownames(plot_detect) <- plot_detect$Freq
  
  plot_list_probe_qc[["percent_segments_by_gene_detection_cutoff"]] <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
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
  
  # Subset to target genes detected in at least x% of the samples
  # Also manually include the negative control probe
  # and any probes that are not in the main module
  neg_probes <- subset(fData(target_data_object), CodeClass == "Negative") %>% .$TargetName %>% unique
  non_main_mod_probes <- subset(fData(target_data_object), Module != main_module) %>% .$TargetName
  target_data_object <- 
    target_data_object[fData(target_data_object)$DetectionRate >= (percent_of_segments/100) | # percent_of_segments is given as a percentage, not a decimal, so we need to divide by 100 first
                         fData(target_data_object)$TargetName %in% union(neg_probes, non_main_mod_probes)
                       , ]
  dim(target_data_object)
  
  # # Retain only detected genes of interest
  # goi <- goi[goi %in% rownames(target_data_object)]
  
  ## ................................................
  ##
  ### Finalization ----
  ##
  ## ................................................
  # Add the main module data object back to the list
  target_data_object_list[[module]] <- target_data_object
  
  # Now subset the segments for the rest of the modules
  for(module_i in names(target_data_object_list)) {
    if(module_i == module) next
    target_data_object_list[[module_i]] <- target_data_object_list[[module_i]][, segments_keep]
  }
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Save the QC-ed NanoStringGeoMxSet object
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-probes.rds"))
# Save the raw probe QC data (fData(data_object))
# We can't write to a table because we have a data-frame-within-a-data-frame
saveRDS(data_object_fData, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-probes-raw.rds"))
# Save the plots
saveRDS(plot_list_probe_qc, paste0(output_dir_rdata, "qc-probes_plot_list.rds"))
# Save the QC table
saveRDS(qc_df, paste0(output_dir_rdata, "qc-probes_table.rds"))
# Save the gene-of-interest detection-rate table
if(exists("goi_df")) saveRDS(goi_df, paste0(output_dir_rdata, "genes-of-interest_detection-rate_table.rds"))
# Save the protein QC table
if(exists("protein_qc_table")) saveRDS(protein_qc_table, paste0(output_dir_rdata, "protein-QC-table.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_qc_probes.RData"))

# Update latest module completed
updateLatestModule(output_dir_rdata, current_module)