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
probes_include_16s <- probes_include_16s %>% str_split(",") %>% unlist()
probes_exclude_16s <- probes_exclude_16s %>% str_split(",") %>% unlist()

# Convert any "NA" in `subset_vars_16s` to "Complete data set"
subset_vars_16s[subset_vars_16s=="NA"] <- "Complete data set"

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
      # Include only the probes given by probes_include_16s (of course including the negative probes)
      target_data_object_16s <- subset(target_data_object_16s, TargetName %in% c(probes_include_16s, neg_probes))
    } else if(str_to_lower(method_16s_probe_selection) == "exclude") {
      # Exclude the probes given by probes_exclude_16s and keep everything else
      target_data_object_16s <- subset(target_data_object_16s, !(TargetName %in% probes_exclude_16s))
    } 
    # Else: don't subset target data object
    
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
    
    # Create a list to hold the excluded levels for each grouping variable/first fixed effect
    # These will be generated in the loop for the differential analysis
    # and then referenced in the visualization loop
    excluded_levels_list <- list()
    
    # Check if the formula table has been provided
    # If formula_table_16s_file is provided, check that it's a valid file
    # If not, skip differential analysis
    valid_formula_table <- TRUE
    if(!flagVariable(formula_table_file_16s)) {
      if(!file.exists(formula_table_file_16s)) {
        warning("The path to the formula table provided does not exist. Differential analysis will not be performed")
        valid_formula_table <- FALSE
      } else {
        # Read in the file and check that it has at least one entry
        message("Checking provided formula table file")
        if(base::grepl("\\.csv$", formula_table_file_16s)) { formula_table_16s <- read.csv(formula_table_file_16s, header = F)}
        else if(base::grepl("\\.tsv$", formula_table_file_16s)) { formula_table_16s <- read.table(formula_table_file_16s, header = F, sep = "\t") }
        else if(base::grepl("\\.xls.*$", formula_table_file_16s)) { formula_table_16s <- read_excel(formula_table_file_16s, col_names = F, na = "NA")}
        else {
          warning("The path to the formula table provided does not exist. Differential analysis will not be performed")
          valid_formula_table <- FALSE
        }
      }
    } else {
      warning("No formula table file was provided. Differential analysis will not be performed")
      valid_formula_table <- FALSE
    }
    # Check that the formula table as read in is valid
    if(exists("formula_table_16s")) if(is.null(formula_table_16s)) {warning("The path provided does not contain a valid formula table. Differential analysis will not be performed"); valid_formula_table <- FALSE}
    
    # If everything is in place to do differential analysis, run it
    if(valid_formula_table) {
      da_res_list <- list()
      
      # If subset variables are not defined, create a dummy subset variable
      subset_vars_16s_orig_flagged <- FALSE # We'll use this in the plotting section, because we need to know if `subset_vars_16s` was originally NULL/empty in the config YAML file
      if(flagVariable(subset_vars_16s)) { subset_vars_16s <- "Complete data set"; subset_vars_16s_orig_flagged <- TRUE }
      
      # Loop over each formula
      # 1 formula / row in formula table
      for(i in 1:nrow(formula_table_16s)) {
        model_number <- i
        da_res_list[[paste0("model_", model_number)]] <- list()
        
        # Extract the formula, whether or not to do all pairwise comparisons, which level to set as the baseline level (if applicable), and any levels of the first fixed effect to exclude
        formula <- formula_table_16s[i,1] %>% as.character
        all_pairwise <- formula_table_16s[i,2]
        baseline_level <- formula_table_16s[i,3] %>% as.character
        excluded_levels <- formula_table_16s[i,4] %>% as.character %>% str_split(";") %>% unlist
        
        # Strip anything before the `~`
        formula <- formula %>% pipe.gsub("^([[:space:]]|.)*~", "~")
        # Extract variables from formula
        formula_vars <- extractVariables(formula) # Input: string
        # Extract first fixed effect from formula
        first_fixed_effect <- extractFirstFixedEffect(as.formula(formula)) # Input: formula
        # Add the dependent variable
        formula_af <- paste0("score ", formula) %>% as.formula
        
        # Now that we have both the first fixed effect (grouping variable) and associated excluded levels,
        # add to `excluded_levels_list` 
        excluded_levels_list[[first_fixed_effect]] <- excluded_levels
        
        # Ensure `Sample` is included
        formula_vars <- c("Sample", formula_vars)
        
        # Add the `Sample` column (created from rownames) to pData, then
        # convert all selected columns (except Sample) to factors
        # and subset to include only necessary variables.
        # First identify which columns in pData_subset are data frames (e.g. LOQ)
        df_cols <- sapply(pData(target_data_object_16s), is.data.frame)
        
        # Create a temporary pData data frame by converting only non-data-frame columns (except `Sample`) to factors
        pData_tmp <- pData(target_data_object_16s) %>%
          tibble::rownames_to_column(var = "Sample") %>% 
          dplyr::mutate(across(!is.data.frame, ~ if (is.numeric(.)) . else as.factor(.)))
        # %>% select(all_of(formula_vars))
        # Add "Complete data set" as variable if it doesn't exist already
        if(("Complete data set" %in% subset_vars_16s) & !("Complete data set" %in% colnames(pData_tmp))) pData_tmp$`Complete data set` <- as.factor("Complete data set")
        
        # If `all_pairwise` is FALSE (i.e., we're doing comparisons against a baseline level)
        # then re-order the factor levels of `first_fixed_effect` in `pData_tmp` so that `baseline_level` is the first
        if(is.logical(all_pairwise)) all_pairwise <- T
        if(!all_pairwise & (baseline_level %in% levels(pData_tmp[[first_fixed_effect]]))) {
          non_baseline_levels <- levels(pData_tmp[[first_fixed_effect]]) %>% .[. != baseline_level]
          new_order <- c(baseline_level, non_baseline_levels)
          pData_tmp[[first_fixed_effect]] <- factor(pData_tmp[[first_fixed_effect]], levels = new_order)
        }
        
        # Merge 16S expression with `pData_tmp` to create the data frame for LMM
        df <- exprs_16s %>% as.data.frame %>% 
          tibble::rownames_to_column(var = "Sample") %>% 
          dplyr::rename(score = ".") %>% # The dependent variable in the formula must be named "score"
          dplyr::left_join(pData_tmp, by = "Sample")
        
        # Loop level 2 (subset variable)
        for(subset_var in subset_vars_16s) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_16s)) subset_vars_16s <- "Complete data set"` above
          if(subset_var == first_fixed_effect) next
          da_res_list[[paste0("model_", model_number)]][[subset_var]] <- list()
          
          if(flagVariable(subset_var) | subset_var == "Complete data set") {
            subset_tag <- "All ROIs"
            subset_var <- "Complete data set"
          } else {
            subset_tag <- subset_var
          }
          
          # Get the levels of the current `subset_var`
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
            # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
            if(!flagVariable(excluded_levels) & !is.na(excluded_levels)) df_final <- df_final %>% dplyr::filter(!(!!as.name(first_fixed_effect) %in% excluded_levels))
            
            # Fit the user-defined linear mixed model
            model <- tryCatch(
              lmerTest::lmer(as.formula(formula_af), data = df_final),
              error = function(e) {
                warning(glue::glue("Error in fitting linear mixed model for model {formula} - {subset_var} - {subset_var_level}: {e$message}"))
                skip_to_next <<- TRUE
              }
            )
            
            if (skip_to_next) next
            message("Successfully calculated differential score")
            
            # If `all_pairwise` is TRUE, extract estimated marginal means (EMMs) for pairwise comparisons
            if(all_pairwise) {
              # Create a dynamic formula for emmeans
              emm_formula <- as.formula(paste0("~ ", first_fixed_effect))
              emm <- emmeans::emmeans(model, emm_formula)  # Replace 'Group' with your fixed effect variable
              # Perform all pairwise comparisons (Tukey-adjusted)
              pairwise_results <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
                as.data.frame()
              
              # Manually create `model_summary` data frame based on `tidy(model)`
              # Colnames of `model_summary` as created by `tidy(model)`: 
              # [1] "effect"           "fixed_effect"     "baseline"         "term"             "estimate"         "std.error"        "statistic"       
              # [8] "df"               "p.value"          "subset_var"       "subset_var_level" "formula"    
              # Colnames of `pairwise_results`:
              # [1] "contrast" "estimate" "SE"       "df"       "t.ratio"  "p.value"
              # So, to turn `pairwise_results` into `model_summary`:
              # add "effect" = "fixed"
              # add "fixed_effect" = `first_fixed_effect`
              # add "baseline" and "term" by splitting "contrast" by "-" and swapping the order
              # on that note, "estimate" in `pairwise_results` is the inverse of "estimate" in `tidy(model)`, so multiply by -1 to get the correct value
              # "std.error" is just "SE", so just rename "SE"
              # like "estimate", "statistic" should be the inverse of "t.ratio" (rename and multiply by -1)
              # "df" is just "df" (although the values are different between the output of `lmer` and `emmeans`)
              # "p.value" is just "p.value" (although, again, the values are different due to multiple testing correction)
              # "subset_var", "subset_var_level", and "formula" need to be added manually
              model_summary <- pairwise_results %>% 
                dplyr::mutate(effect = "fixed", fixed_effect = first_fixed_effect) %>% 
                tidyr::separate(col = "contrast", into = c("baseline", "term"), sep = " - ") %>% 
                dplyr::mutate(estimate = -1 * estimate, statistic = -1 * t.ratio) %>% 
                dplyr::rename(std.error = SE) %>% 
                dplyr::mutate(subset_var = subset_var, subset_var_level = subset_var_level, formula = formula %>% pipe.gsub("~ ", ""))
              model_summary <- model_summary[,c("effect", "fixed_effect", "baseline", "term", "estimate", "std.error",        "statistic", "df", "p.value", "subset_var", "subset_var_level", "formula")]
              
            } else { # Comparisons against a single baseline
              
              # Extract fixed effect results and create `model_summary` data frame
              base_value <- rownames(contrasts(pData_tmp[[first_fixed_effect]]))[1] # Get the baseline value
              model_summary <- tryCatch(
                tidy(model, effects = "fixed") %>%
                  dplyr::filter(term != "(Intercept)") %>% # Filter out the intercept
                  dplyr::mutate(baseline = base_value,
                                subset_var = subset_var,
                                subset_var_level = subset_var_level,
                                fixed_effect = first_fixed_effect,
                                formula = as.character(formula_af)[3]) %>%  # Add cell type column. baseline = base_value, 
                  dplyr::mutate(term = base::gsub(paste0("^", fixed_effect), "", term)) %>% # Clean up the fixed effect name
                  relocate(fixed_effect, .before = term) %>% relocate(baseline, .before = term), # baseline, 
                error = function(e) {
                  warning(glue::glue("Error in extracting fixed effect results for model {formula} - {subset_var} - {subset_var_level}: {e$message}"))
                  skip_to_next <<- TRUE
                }
              )
              
            } # End else comparisons against a single baseline
            if (skip_to_next) next
            message("Successfully extracted fixed effect results")
            
            # Store results
            da_res_list[[paste0("model_", model_number)]][[subset_var]][[subset_var_level]] <- model_summary
            
          } # End loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (LMM model)
          
        } # End loop level 2 (subset variable) << loop level 1 (LMM model)
        
        # # Combine results into a data frame
        # results_df <- bind_rows(da_res_list[[paste0("model_", model_number)]][[method]])
        
      } # End loop level 1 (LMM model)
      
      da_res_df <- bind_rows(rlang::squash(da_res_list)) # `squash` recursively flattens the list
      
    } # End valid formula table check
    
    ## ................................................
    ##
    ### Visualization ----
    ##
    ## ................................................
    if(!is.null(grouping_vars_16s) & (sum(grouping_vars_16s == "") < length(grouping_vars_16s))) {
      
      # Set `subset_vars_16s` to whatever it was in the config YAML file
      if(subset_vars_16s_orig_flagged) subset_vars_16s <- NULL
      # See if we need to do any subsetting prior to graphing
      if(flagVariable(subset_vars_16s) | "Complete data set" %in% subset_vars_16s) {
        # If there are no subset variables/one of them is "Complete data set", we will add a column that will act as a dummy subset variable
        # and change `subset_vars_16s` to be the name of this dummy subset variable
        # This will allow us to use one loop for either case (controls switch 1a or 1b)
        pData(target_data_object_16s)[["Complete data set"]] <- "Complete data set"
        pData(target_data_object_16s)[["Complete data set"]] <- as.factor(pData(target_data_object_16s)[["Complete data set"]])
        subset_vars_16s <- union(subset_vars_16s, "Complete data set") # [is.na(subset_vars_16s)]
        
      } # End control switch 1a (no subset variables) << loop level 1 (model)
      
      # Graph 16S expression levels by group,
      # subsetting if requested
      plot_list <- list()
      
      for(subset_var in subset_vars_16s) {
        plot_list[[subset_var]] <- list()
        
        # Check if it's NA
        if(is.na(subset_var) || subset_var=="NA") {
          subset_var <- "Complete data set"
        }
        
        # Get the levels
        subset_var_levels <- pData(target_data_object_16s)[[subset_var]] %>% unique
        
        # Loop over the levels and subset by each level
        for(subset_var_level in subset_var_levels) {
          plot_list[[subset_var]][[subset_var_level]] <- list()
          
          pData_sub <- pData(target_data_object_16s) %>% dplyr::filter(!!as.name(subset_var) == subset_var_level)
          exprs_16s_sub <- exprs_16s %>% .[names(.) %in% rownames(pData_sub)]
          
          for(grouping_var in grouping_vars_16s) {
            if(grouping_var==subset_var) next
            plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- list()
            
            if(identical(names(exprs_16s_sub), rownames(pData_sub))) {
              plot_df <- data.frame(
                `16S expression` = exprs_16s_sub,
                Group = pData_sub[[grouping_var]]
              )
              colnames(plot_df) <- c("16S expression", grouping_var)
              
              # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
              excluded_levels <- excluded_levels_list[[grouping_var]]
              if(!flagVariable(excluded_levels) & !is.na(excluded_levels)) plot_df <- plot_df %>% dplyr::filter(!(!!as.name(grouping_var) %in% excluded_levels))
              
              # Make sure we have enough colors
              n_colors <- plot_df[[grouping_var]] %>% unique %>% length
              mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
              
              # Plot
              plot_title <- ifelse(subset_var == "Complete data set", "", glue::glue("Subset variable: {subset_var} | level: {subset_var_level}"))
              plot <- plot_df %>%
                ggplot(aes(x = !!as.name(grouping_var), y = `16S expression`, fill = !!as.name(grouping_var))) + 
                geom_boxplot(aes(fill = !!as.name(grouping_var)), 
                             width = 0.3, 
                             lwd = 0.6, 
                             alpha = 0.7, 
                             outlier.shape = NA, 
                             staplewidth = 0.3) + 
                geom_jitter(aes(color = !!as.name(grouping_var)), 
                            width = 0.15, 
                            size = 2, 
                            alpha = 0.9) + 
                scale_fill_manual(values = mycolors) + # , guide = FALSE
                scale_color_manual(values = mycolors) + # , guide = FALSE
                theme(
                  panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.7),
                  axis.line = element_line(linewidth = 0),
                  axis.line.x.top = element_line(linewidth = 0),   # top; only active if we have scales there, e.g. via scale_x_continuous(position = "top")
                  axis.line.y.right = element_line(linewidth = 0), # right; only active if we have scales there
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.position = "right",
                  panel.background = element_rect(fill = "#FFFFFF"),
                  strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
                ) +
                labs(title = plot_title,
                     x = paste0(""),
                     y = paste0("16S expression"))
              
              # If `add_boxplot_pvals_16s` is TRUE, 
              # build p-value data frame and add p-values to plot
              # using `stat_pvalue_manual()`
              if(add_boxplot_pvals_16s) {
                # Check if the current `grouping_var` is in `da_res_df` (the data frame with all the differential LMM results)
                if(grouping_var %in% da_res_df$fixed_effect) {
                  # Create the p-value data frame
                  # p-value data frame: group1, group2, p, y.position, p.label
                  pvals_df <- da_res_df %>%
                    dplyr::filter(fixed_effect == grouping_var & subset_var == !!subset_var & subset_var_level == !!subset_var_level) %>%
                    dplyr::select(baseline, term, p.value) %>%
                    dplyr::rename(group1 = baseline, group2 = term, p = p.value) %>%
                    # Because the LMM or whatever reformats values with special characters
                    # by enclosing them in parentheses ("()"), we need to strip those parentheses
                    # out of the values, otherwise it will cause weird stuff to happen with the graphing
                    dplyr::mutate(group1 = group1 %>% pipe.gsub("^\\(", "") %>% pipe.gsub("\\)$", ""),
                                  group2 = group2 %>% pipe.gsub("^\\(", "") %>% pipe.gsub("\\)$", ""))
                  # If there are no entries with the combination `grouping_var` x `subset_var` x `subset_var_level`, then skip
                  if(nrow(pvals_df) < 1) {message(glue::glue("No entries for combination {grouping_var}, {subset_var}, {subset_var_level}. Skipping")); next}
                  
                  # Calculate y.position: slightly above the highest point per facet
                  tops <- plot_df %>% # `plot_df` should already have only the samples with the correct `subset_var` and `subset_var_level`
                    summarise(y.position = max(`16S expression`, na.rm = TRUE) * 1.05, .groups = "drop")
                  
                  # Add `y.position` to `pvals_df`
                  # The space between brackets should be ~ 10% the range of the points
                  bracket_spacing <- 0.1 * diff(range(plot_df$`16S expression`))
                  highest_bracket <- bracket_spacing * ((pvals_df %>% nrow()) - 1)
                  pvals_df[["y.position"]] <- tops[,1] + seq(from = 0, to = highest_bracket, by = bracket_spacing)
                  # Add `label` column 
                  # pvals_df[["p.label"]] <- sprintf("p = %.2f", pvals_df$p) # %.2f = 2 decimal places; %.2g = 2 significant figures
                  pvals_df <- pvals_df %>% 
                    dplyr::mutate(p.label = dplyr::case_when(
                      p < 0.0001 ~ "****",
                      p < 0.005 ~ "***",
                      p < 0.01 ~ "**",
                      p < 0.05 ~ "*",
                      TRUE ~ sprintf("p = %.2f", p)
                    ))
                  
                  # Add p-values to plot using `ggpubr::stat_pvalue_manual()`
                  plot <- plot + ggpubr::stat_pvalue_manual(
                    pvals_df,
                    y.position = "y.position",
                    label = "p.label",
                    xmin = "group1",
                    xmax = "group2",
                    bracket.size = 0.4,
                    tip.length = 0.01,
                    inherit.aes = FALSE
                  )
                  
                } # End check grouping variable is present in LMM results
              } # End check whether user specified to add p-vals to box plots
              
              # Add to the plot list
              plot_list[[subset_var]][[subset_var_level]][[grouping_var]] <- plot
              # Save plot to disk
              plot
              file_name <- glue::glue("16S_exprs_plot_{subset_var}_{subset_var_level}_{grouping_var}")
              for(file_type in output_plot_file_types) {
                ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 8, height = 6, units = "in")
              }
              
            } else {
              warning("The names in the 16S expression matrix do not match those in the pData. Skipping plots of 16S expression levels by group")
            }
            
          } # End grouping variables `for()` loop
        } # End subset variable levels `for()` loop
      } # End subset variables `for()` loop
      
      # Arrange the plots into a grid
      for(subset_var in names(plot_list)) {
        for(subset_var_level in names(plot_list[[subset_var]])) {
          p_list <- plot_list[[subset_var]][[subset_var_level]]
          
          if(is.null(p_list)) next 
          
          # Remove empty elements from p_list
          p_list <- p_list %>% .[lengths(.) > 0]
          
          n <- length(p_list)
          nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
          plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
          plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("16S expression", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                                                             bottom = grid::textGrob("Grouping variable", gp = grid::gpar(cex = 1.3)),
                                                             top = grid::textGrob(paste0("")))
          
        } # End subset variable level `for()` loop
      } # End subset variable `for()` loop
      
    } # End check for grouping variables
    
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
# Export the LMM results
if(exists("da_res_df")) da_res_df %>% saveRDS(paste0(output_dir_rdata, "16S-analysis_LMM-results.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_16s_analysis.RData"))

# Update latest module completed
updateLatestModule(output_dir_rdata, current_module)