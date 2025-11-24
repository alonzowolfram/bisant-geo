## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
target_data_object_list <- readRDS(cl_args[5])

# Set the probes to be included/excluded
probes_include <- probes_include %>% str_split(",") %>% unlist()
probes_exclude <- probes_exclude %>% str_split(",") %>% unlist()

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## 16S analysis ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(!flagVariable(module_16s) && module_16s %in% names(target_data_object_list)) { # Only run the module if the 16S PKC module is provided & it's in the target_data_object_list
  
  ## ................................................
  ##
  ### Subsetting ----
  ##
  ## ................................................
  # Subset to include only the 16S probes
  target_data_object_16s <- target_data_object_list[[module_16s]] # subset(target_data_object, Module %in% module_16s)
  
  # Stop if there are no 16S probes
  if(!(all(dim(target_data_object_16s) > 0))) {
    warning(glue::glue("There are no probes matching those given in {module_16s}. 16S analysis will not be performed"))
  } else {
    # Get the names of the negative probes
    neg_probes <- fData(target_data_object_16s) %>% dplyr::filter(Negative==TRUE) %>% dplyr::select(TargetName) %>% unlist()
    # Depending on the value given by `method_16s_probe_selection`,
    # exclude or include probes as appropriate
    if(str_to_lower(method_16s_probe_selection) == "include") {
      # Include only the probes given by probes_include (of course including the negative probes)
      target_data_object_16s <- subset(target_data_object_16s, TargetName %in% c(probes_include, neg_probes))
    } else {
      # Exclude the probes given by probes_exclude and keep everything else
      target_data_object_16s <- subset(target_data_object_16s, !(TargetName %in% probes_exclude))
    }
    
    ## ................................................
    ##
    ### Normalization ----
    ##
    ## ................................................
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
    
    # Save to the target data object list
    target_data_object_list[[module_16s]] <- target_data_object_16s
    
    ## ................................................
    ##
    ### 16S score calculation ----
    ##
    ## ................................................
    # The formula for the 16S rRNA score for a given sample: e^((1/n) * sum(ln(x_i)))
    # Where n is the total number of probes,
    # and x_i is the normalized bg-subtracted value for probe i
    bis_mat <- target_data_object_16s@assayData$neg_norm_bgsub %>% .[!(rownames(.) %in% neg_probes),] # Remove the negative probes.
    # Samples in cols, probes in rows.
    # Log transform.
    bis_mat_ln <- log(bis_mat + 1)
    # Sum across probes and divide by the number of probes to get average signal
    n_probes <- nrow(bis_mat)
    mean_signal <- colSums(bis_mat_ln, na.rm = TRUE) / n_probes
    # Exponentiate average signal to get the 16S score
    score_16s <- exp(mean_signal)
    # Add the 16S score to pData for all modules (including 16S)
    for(module in names(target_data_object_list)) {
      if(identical(sampleNames(target_data_object_list[[module]]), names(score_16s)))  {
        pData(target_data_object_list[[module]])$Score16S <- score_16s
      } else {
        warning(glue::glue("The sampleNames for module {module} do not match the names of the 16S score vector. 16S scores will not be added to the pData for this module"))
      }
    }
    
    ## ................................................
    ##
    ### 16S classification ----
    ##
    ## ................................................
    # To determine 16S classification, we'll get the average 16S probe expression
    # for each sample. 
    # Samples will then be categorized as high or low 16S based on (a) user-defined quantile cutoff(s)
    mean_16s_raw <- colMeans(bis_mat)
    if(classification_16s_type=="raw") {exprs_16s <- mean_16s_raw} else {exprs_16s <- score_16s}
    
    percentile_16s_cutoff <- as.numeric(percentile_16s_cutoff)
    for(cutoff in percentile_16s_cutoff) {
      # Determine the classification.
      group_16s <- ifelse(exprs_16s >= quantile(exprs_16s, probs = (cutoff/100), na.rm = TRUE), "16S high", "16S low")
      # Add the 16S scores to pData for all modules (including 16S)
      group_var_name <- paste0("Grouping16S_", cutoff)
      
      for(module in names(target_data_object_list)) {
        pData(target_data_object_list[[module]])[[group_var_name]] <- group_16s
      }
    }
    
    ## ................................................
    ##
    ### Differential 16S levels ----
    ##
    ## ................................................
    if(!flagVariable(lmm_formulae_16s)) {
      model_number <- 1
      da_res_list <- list()
      
      if(flagVariable(subset_vars_16s)) subset_vars_16s <- "Complete dataset"
      
      for(formula in lmm_formulae_16s) {
        da_res_list[[paste0("model_", model_number)]] <- list()
        
        # Strip anything before the `~`.
        formula <- formula %>% regexPipes::gsub("^([[:space:]]|.)*~", "~")
        # Extract variables from formula
        formula_vars <- extractVariables(formula) # Input: string
        # Extract first fixed effect from formula
        first_fixed_effect <- extractFirstFixedEffect(as.formula(formula)) # Input: formula
        # Add the dependent variable
        formula_af <- paste0("score ", formula) %>% as.formula
        
        # Ensure Sample is included
        formula_vars <- c("Sample", formula_vars)
        
        # Add the `Sample` column (created from rownames) to pData, then
        # convert all selected columns (except Sample) to factors
        # and subset to include only necessary variables.
        # First identify which columns in pData_subset are data frames (e.g. LOQ)
        df_cols <- sapply(pData(target_data_object_16s), is.data.frame)
        
        # Convert only non-data-frame columns (except Sample) to factors
        pData_tmp <- pData(target_data_object_16s) %>%
          tibble::rownames_to_column(var = "Sample") %>% 
          dplyr::mutate(across(!is.data.frame, ~ if (is.numeric(.)) . else as.factor(.)))
        # %>% select(all_of(formula_vars))
        # Add "Complete dataset" as variable if it doesn't exist already
        if(subset_vars_16s=="Complete dataset" & !("Complete dataset" %in% colnames(pData_tmp))) pData_tmp$`Complete dataset` <- as.factor("Complete dataset")
        
        # Merge 16S expression with `pData_subset`
        df <- exprs_16s %>% as.data.frame %>% 
          tibble::rownames_to_column(var = "Sample") %>% 
          dplyr::rename(score = ".") %>% # The dependent variable in the formula must be named "score"
          dplyr::left_join(pData_tmp, by = "Sample")
        
        # Loop level 2 (subset variable)
        for(subset_var in subset_vars_16s) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_16s)) subset_vars_16s <- "Complete dataset"` above. 
          da_res_list[[paste0("model_", model_number)]][[subset_var]] <- list()
          
          if(flagVariable(subset_var)) {
            subset_tag <- "All ROIs"
            subset_var <- "Complete dataset"
          } else {
            subset_tag <- subset_var
          }
          
          # Get the levels of the current subset_var.
          subset_var_levels <- pData_tmp[[subset_var]] %>% levels # Previously as.factor %>% needed because it might be a character vector, but now we convert all non-data-frame columns (except Sample) to factor above
          
          # If `subset_var_16s_levels_manual` is set, filter `subset_var_levels` to include only those values
          subset_var_16s_levels_manual_i <- subset_var_16s_levels_manual[[subset_var]]
          if(sum(is.na(subset_var_16s_levels_manual_i)) < length(subset_var_16s_levels_manual_i)) { # At least one `subset_var_level_manual_i` is not NA
            if(sum(subset_var_16s_levels_manual_i %in% subset_var_levels) < 1) {
              warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
            } else { # At least one `subset_var_level_manual_i` is an actual level of the current subset variable
              subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_16s_levels_manual_i]
            }
          }
          
          # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (LMM model)
          for(subset_var_level in subset_var_levels) {
            message("")
            message(glue::glue("Calculating differential abundance using formula {formula} | subset variable {subset_var} | level {subset_var_level} \n"))
            skip_to_next <- FALSE
            
            da_res_list[[paste0("model_", model_number)]][[subset_var]][[subset_var_level]] <- list()
            
            # Get all the samples belonging to the current subset_var_level.
            samples <- pData_tmp %>% 
              dplyr::filter(!!as.name(subset_var)==subset_var_level) %>% 
              dplyr::select(Sample) %>% 
              unlist
            
            # Subset data for this subset variable and subset variable level
            df_final <- df %>%
              dplyr::filter(Sample %in% samples)
            
            # Fit the user-defined linear mixed model
            model <- tryCatch(
              lmerTest::lmer(as.formula(formula_af), data = df_final),
              error = function(e) {
                warning(glue::glue("Error in fitting linear mixed model for model {formula} - {method} - {subset_var} - {subset_var_level} - {cell}: {e$message}"))
                skip_to_next <<- TRUE
              }
            )
            
            if (skip_to_next) next
            message("Successfully calculated differential abundance")
            
            # # Extract estimated marginal means (EMMs) for pairwise comparisons
            # # Create a dynamic formula for emmeans
            # emm_formula <- as.formula(paste0("~ ", first_fixed_effect))
            # emm <- emmeans(model, emm_formula)  # Replace 'Group' with your fixed effect variable
            # # Perform all pairwise comparisons (Tukey-adjusted)
            # pairwise_results <- contrast(emm, method = "pairwise", adjust = "tukey") %>%
            #   as.data.frame() %>%
            #   mutate(CellType = cell)  # Add cell type info
            
            # Extract fixed effect results
            base_value <- rownames(contrasts(pData_tmp[[first_fixed_effect]]))[1] # Get the baseline value.
            model_summary <- tryCatch(
              tidy(model, effects = "fixed") %>%
                dplyr::filter(term != "(Intercept)") %>% # Filter out the intercept.
                dplyr::mutate(baseline = base_value,
                              subset_var = subset_var,
                              subset_var_level = subset_var_level,
                              cell_type = cell, 
                              fixed_effect = first_fixed_effect, 
                              method = method, 
                              formula = as.character(formula_af)[3]) %>%  # Add cell type column. baseline = base_value, 
                dplyr::mutate(term = base::gsub(paste0("^", fixed_effect), "", term)) %>% # Clean up the fixed effect name.
                relocate(fixed_effect, .before = term) %>% relocate(baseline, .before = term), # baseline, 
              error = function(e) {
                warning(glue::glue("Error in extracting fixed effect results for model {formula} - {method} - {subset_var} - {subset_var_level} - {cell}: {e$message}"))
                skip_to_next <<- TRUE
              }
            )
            if (skip_to_next) next
            message("Successfully extracted fixed effect results")
            
            # Store results
            da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]][[subset_var_level]][[cell]] <- model_summary
            
          } # End loop level 3 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
          
        } # End loop level 2 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
        
        # # Combine results into a data frame
        # results_df <- bind_rows(da_res_list[[paste0("model_", model_number)]][[method]])
        
        model_number <- model_number + 1
        
      } # End loop level 1 (LMM model)
      
      da_res_df <- bind_rows(rlang::squash(da_res_list)) # `squash` recursively flattens the list
      
    }
    
    ## ................................................
    ##
    ### Plots of 16S expression by group ----
    ##
    ## ................................................
    if(!is.null(grouping_vars_16s) & (sum(grouping_vars_16s == "") < length(grouping_vars_16s))) {
      # See if we need to do any subsetting prior to graphing
      if(flagVariable(subset_vars_16s)) { # sum(is.na(subset_vars_16s)) == length(subset_vars_16s) || sum(subset_vars_16s=="NA", na.rm = T) == length(subset_vars_16s[!is.na(subset_vars_16s)]) || "NA" %in% subset_vars_16s || NA %in% subset_vars_16s
        # Since there are no subset variables, we will add a column that will act as a dummy subset variable
        # and change subset_vars_16s to be the name of this dummy subset variable
        # This will allow us to use one loop for either case (controls switch 1a or 1b)
        pData(target_data_object_16s)[["Complete data set"]] <- "Dummy level"
        pData(target_data_object_16s)[["Complete data set"]] <- as.factor(pData(target_data_object_16s)[["Complete data set"]])
        subset_vars_16s <- c("Complete data set") # [is.na(subset_vars_16s)]
        
      } # End control switch 1a (no subset variables) << loop level 1 (model)
      
      # Graph 16S expression levels by group,
      # subsetting if requested
      plot_list <- list()
      anova_list <- list()
      
      for(subset_var in subset_vars_16s) {
        plot_list[[subset_var]] <- list()
        anova_list[[subset_var]] <- list()
        
        # Check if it's NA
        if(is.na(subset_var) || subset_var=="NA") {
          subset_var <- "Complete data set"
        }
        
        # Get the levels
        subset_var_levels <- pData(target_data_object_16s)[[subset_var]] %>% unique
        
        # Loop over the levels and subset by each level
        for(subset_var_level in subset_var_levels) {
          plot_list[[subset_var]][[subset_var_level]] <- list()
          anova_list[[subset_var]][[subset_var_level]] <- list()
          
          pData_sub <- pData(target_data_object_16s) %>% dplyr::filter(!!as.name(subset_var) == subset_var_level)
          mean_16s_sub <- mean_16s_raw %>% .[names(.) %in% rownames(pData_sub)]
          
          for(grouping_var in grouping_vars_16s) {
            plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            anova_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            
            if(identical(names(mean_16s_sub), rownames(pData_sub))) {
              dat <- data.frame(
                `16S expression` = mean_16s_sub,
                Group = pData_sub[[grouping_var]]
              )
              colnames(dat) <- c("16S expression", grouping_var)
              
              # Check that we have enough levels to perform ANOVA
              n_grouping_var_levels <- dat$Group %>% unique %>% length
              if(n_grouping_var_levels < 2) {
                warning(glue::glue("Grouping variable {grouping_var} has fewer than 2 levels. Skipping differential analysis"))
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- NA
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- NA
              } else {
                # Perform ANOVA.
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- anova_res
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- tukey_res
              }
              
              # Make sure we have enough colors
              n_colors <- dat[[grouping_var]] %>% unique %>% length
              mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
              
              # Plot
              plot <- dat %>%
                ggplot(aes(x = !!as.name(grouping_var), y = `16S expression`, fill = !!as.name(grouping_var))) + 
                geom_boxplot(aes(fill = !!as.name(grouping_var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
                geom_jitter(aes(color = !!as.name(grouping_var)), width = 0.15, size = 2, alpha = 0.9) + 
                #facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
                scale_fill_manual(values = mycolors) + # , guide = FALSE
                scale_color_manual(values = mycolors) + # , guide = FALSE
                theme_bw() +
                theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "right",
                  strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
                ) +
                labs(title = glue::glue("Subset variable: {subset_var} | level: {subset_var_level}"),
                     x = paste0(""),
                     y = paste0("16s expression"))
              # Add to the plot list
              plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- plot
              # Save plot to disk
              plot
              file_name <- glue::glue("16S_plot_{subset_var}_{subset_var_level}_{grouping_var}")
              for(file_type in output_plot_file_types) {
                ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
              }
              
            } else {
              warning("The names in the 16S expression matrix do not match those in the pData! Skipping plots of 16S expression levels by group")
            }
            
          } # End grouping variables for loop
          
        } # End subset variable levels for loop
        
      } # End subset_vars_16s for loop
      
      # Arrange the plots into a grid
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
          
        }
      }
      
    }
    
    ## ................................................
    ##
    ### 16S raw expression levels by group ----
    ##
    ## ................................................
    if(!is.null(grouping_vars_16s) & (sum(grouping_vars_16s == "") < length(grouping_vars_16s))) {
      # See if we need to do any subsetting prior to graphing
      if(flagVariable(subset_vars_16s)) { # sum(is.na(subset_vars_16s)) == length(subset_vars_16s) || sum(subset_vars_16s=="NA", na.rm = T) == length(subset_vars_16s[!is.na(subset_vars_16s)]) || "NA" %in% subset_vars_16s || NA %in% subset_vars_16s
        # Since there are no subset variables, we will add a column that will act as a dummy subset variable
        # and change subset_vars_16s to be the name of this dummy subset variable
        # This will allow us to use one loop for either case (controls switch 1a or 1b)
        pData(target_data_object_16s)[["Complete data set"]] <- "Dummy level"
        pData(target_data_object_16s)[["Complete data set"]] <- as.factor(pData(target_data_object_16s)[["Complete data set"]])
        subset_vars_16s <- c("Complete data set") # [is.na(subset_vars_16s)]
        
      } # End control switch 1a (no subset variables) << loop level 1 (model)
      
      # Graph 16S expression levels by group,
      # subsetting if requested
      plot_list <- list()
      anova_list <- list()
      
      for(subset_var in subset_vars_16s) {
        plot_list[[subset_var]] <- list()
        anova_list[[subset_var]] <- list()
        
        # Check if it's NA
        if(is.na(subset_var) || subset_var=="NA") {
          subset_var <- "Complete data set"
        }
        
        # Get the levels
        subset_var_levels <- pData(target_data_object_16s)[[subset_var]] %>% unique
        
        # Loop over the levels and subset by each level
        for(subset_var_level in subset_var_levels) {
          plot_list[[subset_var]][[subset_var_level]] <- list()
          anova_list[[subset_var]][[subset_var_level]] <- list()
          
          pData_sub <- pData(target_data_object_16s) %>% dplyr::filter(!!as.name(subset_var) == subset_var_level)
          mean_16s_sub <- mean_16s_raw %>% .[names(.) %in% rownames(pData_sub)]
          
          for(grouping_var in grouping_vars_16s) {
            plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            anova_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            
            if(identical(names(mean_16s_sub), rownames(pData_sub))) {
              dat <- data.frame(
                `16S expression` = mean_16s_sub,
                Group = pData_sub[[grouping_var]]
              )
              colnames(dat) <- c("16S expression", grouping_var)
              
              # Check that we have enough levels to perform ANOVA
              n_grouping_var_levels <- dat$Group %>% unique %>% length
              if(n_grouping_var_levels < 2) {
                warning(glue::glue("Grouping variable {grouping_var} has fewer than 2 levels. Skipping ANOVA"))
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- NA
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- NA
              } else {
                # Perform ANOVA.
                anova_res <- aov(`16S expression` ~ get(grouping_var), data = dat)
                tukey_res <- TukeyHSD(anova_res)
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["ANOVA"]] <- anova_res
                anova_list[[subset_var]][[subset_var_level]][[grouping_var]][["TukeyHSD"]] <- tukey_res
              }
              
              # Make sure we have enough colors
              n_colors <- dat[[grouping_var]] %>% unique %>% length
              mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
              
              # Plot
              plot <- dat %>%
                ggplot(aes(x = !!as.name(grouping_var), y = `16S expression`, fill = !!as.name(grouping_var))) + 
                geom_boxplot(aes(fill = !!as.name(grouping_var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
                geom_jitter(aes(color = !!as.name(grouping_var)), width = 0.15, size = 2, alpha = 0.9) + 
                #facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
                scale_fill_manual(values = mycolors) + # , guide = FALSE
                scale_color_manual(values = mycolors) + # , guide = FALSE
                theme_bw() +
                theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "right",
                  strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
                ) +
                labs(title = glue::glue("Subset variable: {subset_var} | level: {subset_var_level}"),
                     x = paste0(""),
                     y = paste0("16s expression"))
              # Add to the plot list
              plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- plot
              # Save plot to disk
              plot
              file_name <- glue::glue("16S_plot_{subset_var}_{subset_var_level}_{grouping_var}")
              for(file_type in output_plot_file_types) {
                ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
              }
              
            } else {
              warning("The names in the 16S expression matrix do not match those in the pData! Skipping plots of 16S expression levels by group")
            }
            
          } # End grouping variables for loop
          
        } # End subset variable levels for loop
        
      } # End subset_vars_16s for loop
      
      # Arrange the plots into a grid
      for(subset_var in names(plot_list)) {
        for(subset_var_level in names(plot_list[[subset_var]])) {
          p_list <- plot_list[[subset_var]][[subset_var_level]]
          
          if(is.null(p_list)) next 
          
          n <- length(p_list)
          nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1
          plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
          plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("16S expression", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                                                             bottom = grid::textGrob("Grouping variable", gp = grid::gpar(cex = 1.3)),
                                                             top = grid::textGrob(paste0("")))
          
        }
      }
      
    }
    
  } # End check for 16S probes
} # End check for 16S module

# Check if there's an additional column added if there are no subset variables
# If there is, remove it
# (Otherwise, this will mess things up downstream)
for(module in names(target_data_object_list)) {
  if("Complete data set" %in% colnames(pData(target_data_object_list[[module]]))) pData(target_data_object_list[[module]]) <- pData(target_data_object_list[[module]]) %>% dplyr::select(-`Complete data set`)
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export NanoStringGeoMxSet as RDS file
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_16S-analysis.rds"))
# Export graphs
if(exists("plot_list")) saveRDS(plot_list, paste0(output_dir_rdata, "16S-analysis_raw-plots-list.rds"))
# Export ANOVA results
if(exists("anova_list")) saveRDS(anova_list, paste0(output_dir_rdata, "16S-analysis_ANOVA-res-list.rds"))

# Update latest module completed
updateLatestModule(output_dir_rdata, current_module)