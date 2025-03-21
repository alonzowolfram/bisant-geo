## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object.
target_data_object_list <- readRDS(cl_args[5])
# We'll only need the 16S module for this one.
if(!flagVariable(module_16s)) {
  target_data_object <- target_data_object_list[[module_16s]]
} else {
  target_data_object <- target_data_object_list[[main_module]]
}

# Set the probes to be included/excluded. 
probes_include <- probes_include %>% str_split(",") %>% unlist()
probes_exclude <- probes_exclude %>% str_split(",") %>% unlist()

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## 16S analysis ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(!flagVariable(module_16s)) { # Only run the module if the 16S PKC module is provided & it's in the target_data_object_list.
  
  ## ................................................
  ##
  ### Subsetting ----
  ##
  ## ................................................
  # Subset to include only the 16S probes.
  target_data_object_16s <- subset(target_data_object, Module %in% module_16s)
  
  # Stop if there are no 16S probes.
  if(!(all(dim(target_data_object_16s) > 0))) {
    warning(paste0("There are no probes matching those given in", module_16s, ". 16S analysis will not be performed."))
  } else {
    # Get the names of the negative probes.
    neg_probes <- fData(target_data_object_16s) %>% dplyr::filter(Negative==TRUE) %>% dplyr::select(TargetName) %>% unlist()
    # Depending on the value given by `method_16s_probe_selection`,
    # exclude or include probes as appropriate.
    if(str_to_lower(method_16s_probe_selection) == "include") {
      # Include only the probes given by probes_include (of course including the negative probes).
      target_data_object_16s <- subset(target_data_object_16s, TargetName %in% c(probes_include, neg_probes))
    } else {
      # Exclude the probes given by probes_exclude and keep everything else.
      target_data_object_16s <- subset(target_data_object_16s, !(TargetName %in% probes_exclude))
    }
    
    ## ................................................
    ##
    ### Normalization ----
    ##
    ## ................................................
    
    # # Subtract background to account for signal:noise discrepancies
    # target_data_object_16s <- normalize(target_data_object_16s,
    #                                  norm_method = "subtractBackground",
    #                                  fromElt = "exprs",
    #                                  toElt = "bgsub"
    # )
    # Already done in normalization module?
    
    # Calculate background normalization factor
    pData(target_data_object_16s)$neg_normFactor <-
      pData(target_data_object_16s)[[paste0("NegGeoMean_", module_16s)]] /
      ngeoMean(pData(target_data_object_16s)[[paste0("NegGeoMean_", module_16s)]])
    
    # Normalize to background
    assayDataElement(target_data_object_16s, "neg_norm_bgsub") <-
      sweep(assayDataElement(target_data_object_16s, "bg_sub"), 2L,
            pData(target_data_object_16s)$neg_normFactor,
            FUN = "/"
      )
    
    ## ................................................
    ##
    ### 16S score ----
    ##
    ## ................................................
    # The formula for the 16S rRNA score for a given sample: e^((1/n) * sum(ln(x_i)))
    # Where n is the total number of probes,
    # and x_i is the normalized bg-subtracted value for probe i.
    bis_mat <- target_data_object_16s@assayData$neg_norm_bgsub %>% .[!(rownames(.) %in% neg_probes),] # Remove the negative probes.
    # Samples in cols, probes in rows.
    # Log transform.
    bis_mat_ln <- log(bis_mat + 1)
    # Sum across probes and divide by the number of probes to get average signal.
    n_probes <- nrow(bis_mat)
    mean_signal <- colSums(bis_mat_ln, na.rm = TRUE) / n_probes
    # Exponentiate average signal to get the 16S score.
    score_16s <- exp(mean_signal)
    # Add the 16S score to the pData.
    if(identical(sampleNames(target_data_object_16s), names(score_16s)))  {
      pData(target_data_object_16s)$Score16S <- score_16s
      pData(target_data_object_list[[module_16s]])$Score16S <- score_16s
    }
    
    # Add the 16S scores to the other pDatas.
    for(module in names(target_data_object_list)) {
      if(module == module_16s) next
      
      pData(target_data_object_list[[module]])$Score16S <- score_16s
    }
    
    ## ................................................
    ##
    ### 16S classification ----
    ##
    ## ................................................
    # To determine 16S classification, we'll get the average 16S probe expression
    # for each sample. 
    # Samples will then be categorized as high or low 16S based on (a) user-defined quantile cutoff(s).
    mean_16s <- colMeans(bis_mat)
    
    percentile_16s_cutoff <- as.numeric(percentile_16s_cutoff)
    for(cutoff in percentile_16s_cutoff) {
      # Determine the classification.
      group_16s <- ifelse(mean_16s >= quantile(mean_16s, probs = (cutoff/100), na.rm = TRUE), "16S high", "16S low")
      # Add to the metadata of the original target_data_object.
      group_var_name <- paste0("Grouping16S_", cutoff)
      pData(target_data_object)[[group_var_name]] <- group_16s
    }
    
    ## ................................................
    ##
    ### 16S expression levels by group ----
    ##
    ## ................................................
    if(!is.null(exprs_16s_grouping_vars) & (sum(exprs_16s_grouping_vars == "") < length(exprs_16s_grouping_vars))) {
      # # Add a section header.
      # pptx <- pptx %>%
      #   officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
      #   officer::ph_with(value = paste0("16S analysis"),
      #                    location = ph_location_label(ph_label = "Title 1"))
      
      # See if we need to do any subsetting prior to graphing.
      if(sum(is.na(exprs_16s_subset_vars)) == length(exprs_16s_subset_vars) || sum(exprs_16s_subset_vars=="NA", na.rm = T) == length(exprs_16s_subset_vars[!is.na(exprs_16s_subset_vars)]) || "NA" %in% exprs_16s_subset_vars || NA %in% exprs_16s_subset_vars) {
        # Since there are no subset variables, we will add a column that will act as a dummy subset variable
        # and change exprs_16s_subset_vars to be the name of this dummy subset variable.
        # This will allow us to use one loop for either case (controls switch 1a or 1b).
        pData(target_data_object)[["Complete data set"]] <- "Dummy level"
        pData(target_data_object)[["Complete data set"]] <- as.factor(pData(target_data_object)[["Complete data set"]])
        exprs_16s_subset_vars[is.na(exprs_16s_subset_vars)] <- c("Complete data set")
        
      } # End control switch 1a (no subset variables) << loop level 1 (model).
      
      # Graph 16S expression levels by group,
      # subsetting if requested.
      plot_list <- list()
      anova_list <- list()
      
      for(subset_var in exprs_16s_subset_vars) {
        plot_list[[subset_var]] <- list()
        anova_list[[subset_var]] <- list()
        
        # Check if it's NA.
        if(is.na(subset_var) || subset_var=="NA") {
          subset_var <- "Complete data set"
        }
        
        # Get the levels. 
        subset_var_levels <- pData(target_data_object)[[subset_var]] %>% unique
        
        # Loop over the levels and subset by each level.
        for(subset_var_level in subset_var_levels) {
          plot_list[[subset_var]][[subset_var_level]] <- list()
          anova_list[[subset_var]][[subset_var_level]] <- list()
          
          pData_sub <- pData(target_data_object) %>% dplyr::filter(!!as.name(subset_var) == subset_var_level)
          mean_16s_sub <- mean_16s %>% .[names(.) %in% rownames(pData_sub)]
          
          for(grouping_var in exprs_16s_grouping_vars) {
            plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            anova_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            
            if(identical(names(mean_16s_sub), rownames(pData_sub))) {
              dat <- data.frame(
                `16S expression` = mean_16s_sub,
                Group = pData_sub[[grouping_var]]
              )
              colnames(dat) <- c("16S expression", grouping_var)
              
              # Check that we have enough levels to perform ANOVA.
              n_grouping_var_levels <- dat$Group %>% unique %>% length
              if(n_grouping_var_levels < 2) {
                warning("Grouping variable ", grouping_var, " has fewer than 2 levels. Skipping ANOVA.")
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- NA
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- NA
              } else {
                # Perform ANOVA.
                anova_res <- aov(`16S expression` ~ get(grouping_var), data = dat)
                tukey_res <- TukeyHSD(anova_res)
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- anova_res
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- tukey_res
              }
              
              # Make sure we have enough colors.
              n_colors <- dat[[grouping_var]] %>% unique %>% length
              mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
              
              # Plot.
              plot <- dat %>%
                ggplot(aes(x = !!as.name(grouping_var), y = `16S expression`, fill = !!as.name(grouping_var))) + 
                geom_boxplot() + 
                #facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
                scale_fill_manual(values = mycolors) + # , guide = FALSE
                scale_color_manual(values = mycolors) + # , guide = FALSE
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank()) +
                labs(title = paste0("Subset variable: ", subset_var, " | level: ", subset_var_level),
                     x = paste0(""),
                     y = paste0(""))
              plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- plot
              
            } else {
              warning("The names in the 16S expression matrix do not match those in the pData! Skipping plots of 16S expression levels by group.")
            }
            
          } # End grouping variables for loop.
          
        } # End subset variable levels for loop.
        
      } # End exprs_16s_subset_vars for loop.
      
      # # Graphing parameters.
      # plot_width <- 20
      # plot_height <- 15
      # units <- "in"
      # res <- 300
      # scaling_factor <- 1
      # res_scaling_factor <- 1
      # Arrange the plots into a grid.
      for(subset_var in names(plot_list)) {
        for(subset_var_level in names(plot_list[[subset_var]])) {
          p_list <- plot_list[[subset_var]][[subset_var_level]]
          
          if(is.null(p_list)) next 
          
          n <- length(p_list)
          nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
          plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
          plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("16S expression", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                                                             bottom = grid::textGrob("Grouping variable", gp = grid::gpar(cex = 1.3)),
                                                             top = grid::textGrob(paste0("")))
          
          
          # # Save to EPS and PNG and then ...
          # eps_path <- paste0(output_dir_pubs, "16S_exprs_by-group.eps") #paste0(output_dir_pubs, "")
          # png_path <- paste0(output_dir_pubs, "16S_exprs_by-group.png") #paste0(output_dir_pubs, "")
          # plot <- plot_grid# %>% ggpubr::as_ggplot()
          # saveEPS(plot, eps_path, width = (plot_width * scaling_factor), height = (plot_height * scaling_factor))
          # savePNG(plot, png_path, width = (plot_width * scaling_factor), height = (plot_height * scaling_factor), units = units, res = (res * res_scaling_factor))
          
          # # ... add to the PowerPoint.
          # pptx <- pptx %>%
          #   officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
          #   officer::ph_with(value = paste0("16S analysis"),
          #                    location = ph_location_label(ph_label = "Title 1")) %>%
          #   officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
          #                    location = ph_location_label(ph_label = "Content Placeholder 2"),
          #                    use_loc_size = FALSE)
          
        }
      }
      
    }
    
    # Add the 16S classifications to the other pDatas.
    for(module in names(target_data_object_list)) {
      if(module == module_16s) next
      
      # Get the columns that the 16s pData has that the other one doesn't.
      missing_cols <- setdiff(colnames(pData(target_data_object)), colnames(pData(target_data_object_list[[module]])))
      # Add them.
      pData(target_data_object_list[[module]]) <- cbind(pData(target_data_object_list[[module]]), pData(target_data_object)[,missing_cols])
    }
    
  } # End check for 16S probes.
} # End check for 16S module.

# Check if there's an additional column added if there are no subset variables.
# If there is, remove it. 
# (Otherwise, this will mess things up downstream.)
if("Complete data set" %in% colnames(pData(target_data_object))) pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-`Complete data set`)

# Now save the 16S data object back to the list.
if(exists("target_data_object_16s") && all(dim(target_data_object_16s) > 0)) target_data_object_list[[module_16s]] <- target_data_object

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export NanoStringGeoMxSet as RDS file.
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_16S-analysis.rds"))
# Export graphs.
if(exists("plot_list")) saveRDS(plot_list, paste0(output_dir_rdata, "16S-analysis_raw-plots-list.rds"))
# Export ANOVA results.
if(exists("anova_list")) saveRDS(anova_list, paste0(output_dir_rdata, "16S-analysis_ANOVA-res-list.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)