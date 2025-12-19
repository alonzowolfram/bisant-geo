## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
target_data_object_list <- readRDS(cl_args[5])

# Convert any "NA" in `subset_vars_tcr` to "Complete data set"
subset_vars_tcr[subset_vars_tcr=="NA"] <- "Complete data set"

# Functions
# Gini coefficient
calculate_gini_coefficient <- function(x) {
  weights <- rep(1 / length(x), length(x))
  
  x <- x[order(x)]
  p <- cumsum(weights)
  n <- length(x)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu / nu[n]
  
  gini_coeff <- sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
  
  return(gini_coeff)
}
# Shannon diversity
calculate_shannon_diversity <- function(x) {
  total_counts <- apply(x, 1, sum)
  x <- sweep(x, 1, total_counts, "/")
  x <- -x * log(x, exp(1))
  H <- apply(x, 1, sum, na.rm = TRUE)
  
  return(H)
}

# Names of the metrics as stored in pData, for when we loop over them later
metric_names <- c("Gini", "ShannonH", "Simpson", "InvSimpson")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## TCR analysis ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if(!flagVariable(module_tcr) && module_tcr %in% names(target_data_object_list)) { # Only run the module if the TCR module is provided and the combined module is in the target data object list
  # Get the TCR object
  target_data_object <- target_data_object_list[[module_tcr]]
  
  # Get TCR probes
  tcr_probes <- fData(target_data_object)$TargetName[base::grepl("TR[A/B/D/G][C/J/V]", fData(target_data_object)$TargetName)]
  if(length(tcr_probes) < 1) { # No TCR probes
    
    warning("There are no TCR probes in the QC-ed data set. No analysis will be performed")
    plot_list <- list()
    plot_grid_list <- list()
    anova_list <- list()
    
  } else { # There are TCR probes
    
    ### ................................................
    ###
    ### TCR diversity and distribution calculations ----
    ###
    ### ................................................
    
    # Calculate distribution (Gini coefficient)
    gini <- apply(assayDataElement(target_data_object, elt = normalization_tcr)[tcr_probes, ], 2, calculate_gini_coefficient)
    gini[!is.finite(gini)] <- NA
    pData(target_data_object)$Gini <- gini
    
    # Calculate diversity (Shannon, Simpson, and inverse Simpson)
    shannon_h <- vegan::diversity(t(assayData(target_data_object)[[normalization_tcr]][tcr_probes, ]), index = "shannon", MARGIN = 1)
    simpson <- vegan::diversity(t(assayData(target_data_object)[[normalization_tcr]][tcr_probes, ]), index = "simpson", MARGIN = 1)
    invsimpson <- vegan::diversity(t(assayData(target_data_object)[[normalization_tcr]][tcr_probes, ]), index = "invsimpson", MARGIN = 1)
    shannon_h[!is.finite(shannon_h)] <- NA
    simpson[!is.finite(simpson)] <- NA
    invsimpson[!is.finite(invsimpson)] <- NA
    pData(target_data_object)$ShannonH <- shannon_h
    pData(target_data_object)$Simpson <- simpson
    pData(target_data_object)$InvSimpson <- invsimpson
    
    ### ................................................
    ###
    ### Differential analysis of diversity and distribution measures ----
    ###
    ### ................................................
    
    # Create a list to hold the excluded levels for each grouping variable/first fixed effect
    # These will be generated in the loop for the differential analysis
    # and then referenced in the visualization loop
    excluded_levels_list <- list()
    
    # Check if the formula table has been provided
    # If formula_table_tcr_file is provided, check that it's a valid file
    # If not, skip differential analysis
    valid_formula_table <- TRUE
    if(!flagVariable(formula_table_file_tcr)) { 
      if(!file.exists(formula_table_file_tcr)) {
        warning("The path to the formula table provided does not exist. Differential analysis will not be performed")
        valid_formula_table <- FALSE
      } else {
        # Read in the file and check that it has at least one entry
        message("Checking provided formula table file")
        if(base::grepl("\\.csv$", formula_table_file_tcr)) { formula_table_tcr <- read.csv(formula_table_file_tcr, header = F)}
        else if(base::grepl("\\.tsv$", formula_table_file_tcr)) { formula_table_tcr <- read.table(formula_table_file_tcr, header = F, sep = "\t") }
        else if(base::grepl("\\.xls.*$", formula_table_file_tcr)) { formula_table_tcr <- read_excel(formula_table_file_tcr, col_names = F, na = "NA")} 
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
    if(exists("formula_table_tcr")) if(is.null(formula_table_tcr)) {warning("The path provided does not contain a valid formula table. Differential analysis will not be performed"); valid_formula_table <- FALSE}
    
    # If everything is in place to do differential analysis, run it
    if(valid_formula_table) {
      da_res_list <- list()
      
      # If subset variables are not defined, create a dummy subset variable
      subset_vars_tcr_orig_flagged <- FALSE # We'll use this in the plotting section, because we need to know if `subset_vars_tcr` was originally NULL/empty in the config YAML file
      if(flagVariable(subset_vars_tcr)) { subset_vars_tcr <- "Complete data set"; subset_vars_tcr_orig_flagged <- TRUE }
      
      # Loop over each formula
      # 1 formula / row in formula table
      for(i in 1:nrow(formula_table_tcr)) {
        model_number <- i
        da_res_list[[paste0("model_", model_number)]] <- list()
        
        # Extract the formula, whether or not to do all pairwise comparisons, which level to set as the baseline level (if applicable), and any levels of the first fixed effect to exclude
        formula <- formula_table_tcr[i,1] %>% as.character
        all_pairwise <- formula_table_tcr[i,2]
        baseline_level <- formula_table_tcr[i,3] %>% as.character
        excluded_levels <- formula_table_tcr[i,4] %>% as.character %>% str_split(";") %>% unlist
        
        # Strip anything before the `~`
        formula <- formula %>% regexPipes::gsub("^([[:space:]]|.)*~", "~")
        # Extract variables from formula
        formula_vars <- extractVariables(formula) # Input: string
        # Ensure Sample is included
        formula_vars <- c("Sample", formula_vars)
        # Extract first fixed effect from formula
        first_fixed_effect <- extractFirstFixedEffect(as.formula(formula)) # Input: formula
        
        # Now that we have both the first fixed effect (grouping variable) and associated excluded levels,
        # add to `excluded_levels_list` 
        excluded_levels_list[[first_fixed_effect]] <- excluded_levels
        
        # 2025/11/26: for 16S analysis, we would add the dependent variable to create `formula_af` 
        # at this point in the workflow
        # However, because we have 4 different dependent variables this time (Gini, Shannon, etc.)
        # we will create a loop and update `formula_af` within the loop, further down the line
        # (once we're inside the subsetting variable loops)
    
        # Add the `Sample` column (created from rownames) to pData, then
        # convert all selected columns (except Sample) to factors
        # and subset to include only necessary variables.
        # First identify which columns in pData_subset are data frames (e.g. LOQ)
        df_cols <- sapply(pData(target_data_object), is.data.frame)
        
        # Create a temporary pData data frame by converting only non-data-frame columns (except Sample) to factors
        pData_tmp <- pData(target_data_object) %>%
          tibble::rownames_to_column(var = "Sample") %>% 
          dplyr::mutate(across(!is.data.frame, ~ if (is.numeric(.)) . else as.factor(.)))
        # %>% select(all_of(formula_vars))
        # Add "Complete data set" as variable if it doesn't exist already
        if(("Complete data set" %in% subset_vars_tcr) & !("Complete data set" %in% colnames(pData_tmp))) pData_tmp$`Complete data set` <- as.factor("Complete data set")
        
        # If `all_pairwise` is FALSE (i.e., we're doing comparisons against a baseline level)
        # then re-order the factor levels of `first_fixed_effect` in `pData_tmp` so that `baseline_level` is the first
        if(is.logical(all_pairwise)) all_pairwise <- T
        if(!all_pairwise & (baseline_level %in% levels(pData_tmp[[first_fixed_effect]]))) {
          non_baseline_levels <- levels(pData_tmp[[first_fixed_effect]]) %>% .[. != baseline_level]
          new_order <- c(baseline_level, non_baseline_levels)
          pData_tmp[[first_fixed_effect]] <- factor(pData_tmp[[first_fixed_effect]], levels = new_order)
        }
        
        # Create the data frame for LMM from pData
        # Unlike with 16S analysis, the dependent variables are already in the pData
        # so we don't need to left-join separate data frames
        df <- pData(target_data_object) %>% tibble::rownames_to_column(var = "Sample") %>% as.data.frame
        
        # Loop level 2 (subset variable)
        for(subset_var in subset_vars_tcr) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_tcr)) subset_vars_tcr <- "Complete data set"` above
          if(subset_var == first_fixed_effect) next
          da_res_list[[paste0("model_", model_number)]][[subset_var]] <- list()
          
          if(flagVariable(subset_var) | subset_var == "Complete data set") {
            subset_tag <- "All ROIs"
            subset_var <- "Complete data set"
          } else {
            subset_tag <- subset_var
          }
          
          # Get the levels of the current `subset_var`
          subset_var_levels <- pData_tmp[[subset_var]] %>% levels # Previously as.factor %>% needed because it might be a character vector, but now we convert all non-data-frame columns (except `Sample`) to factor above
          
          # If `subset_var_tcr_levels_manual` is set, filter `subset_var_levels` to include only those values
          subset_var_tcr_levels_manual_i <- subset_var_tcr_levels_manual[[subset_var]]
          if(sum(is.na(subset_var_tcr_levels_manual_i)) < length(subset_var_tcr_levels_manual_i)) { # At least one `subset_var_level_manual_i` is not NA
            if(sum(subset_var_tcr_levels_manual_i %in% subset_var_levels) < 1) {
              warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
            } else { # At least one `subset_var_level_manual_i` is an actual level of the current subset variable
              subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_tcr_levels_manual_i]
            }
          }
          
          # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (LMM model)
          for(subset_var_level in subset_var_levels) {
            message("")
            message(glue::glue("Calculating differential abundance using formula {formula} | subset variable {subset_var} | level {subset_var_level} \n"))
            skip_to_next <- FALSE
            
            da_res_list[[paste0("model_", model_number)]][[subset_var]][[subset_var_level]] <- list()
            
            # Get all the samples belonging to the current `subset_var_level`
            samples <- pData_tmp %>%
              dplyr::filter(!!as.name(subset_var)==subset_var_level) %>%
              dplyr::select(Sample) %>%
              unlist
            
            # Subset data for this subset variable and subset variable level
            df_final <- df %>%
              dplyr::filter(Sample %in% samples)
            # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
            if(!flagVariable(excluded_levels)) df_final <- df_final %>% dplyr::filter(!(!!as.name(first_fixed_effect) %in% excluded_levels))
            
            # Fit the user-defined linear mixed model
            # Loop over the different diversity/distribution metrics
            for(metric in metric_names) {
              # Add the dependent variable
              formula_af <- paste(metric, formula) %>% as.formula
              
              model <- tryCatch(
                expr = lmerTest::lmer(as.formula(formula_af), data = df_final),
                error = function(e) {
                  warning(glue::glue("Error in fitting linear mixed model for model {formula} - {subset_var} - {subset_var_level}: {e$message}"))
                  skip_to_next <<- TRUE
                }
              )
              
              if (skip_to_next) next
              message(glue::glue("Successfully calculated differential score for metric {metric}"))
              
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
                  dplyr::mutate(subset_var = subset_var, 
                                subset_var_level = subset_var_level, 
                                formula = formula %>% regexPipes::gsub("~ ", ""),
                                metric = metric)
                model_summary <- model_summary[,c("metric", "effect", "fixed_effect", "baseline", "term", "estimate", "std.error",        "statistic", "df", "p.value", "subset_var", "subset_var_level", "formula")]
                
              } else { # Comparisons against a single baseline
                
                # Extract fixed effect results and create `model_summary` data frame
                base_value <- rownames(contrasts(pData_tmp[[first_fixed_effect]]))[1] # Get the baseline value
                model_summary <- tryCatch(
                  expr = tidy(model, effects = "fixed") %>%
                    dplyr::filter(term != "(Intercept)") %>% # Filter out the intercept
                    dplyr::mutate(metric = metric,
                                  baseline = base_value,
                                  subset_var = subset_var,
                                  subset_var_level = subset_var_level,
                                  fixed_effect = first_fixed_effect,
                                  formula = as.character(formula_af)[3]) %>%  # Add cell type column. baseline = base_value, 
                    dplyr::mutate(term = base::gsub(paste0("^", fixed_effect), "", term)) %>% # Clean up the fixed effect name
                    dplyr::relocate(fixed_effect, .before = term) %>% 
                    dplyr::relocate(baseline, .before = term) %>%
                    dplyr::relocate(metric, .before = effect), # baseline, 
                  error = function(e) {
                    warning(glue::glue("Error in extracting fixed effect results for model {formula} - {subset_var} - {subset_var_level}: {e$message}"))
                    skip_to_next <<- TRUE
                  }
                ) # End tryCatch() for model summary
                
              } # End else comparisons against a single baseline
              if (skip_to_next) next
              message("Successfully extracted fixed effect results")
              
              # Store results
              da_res_list[[paste0("model_", model_number)]][[subset_var]][[subset_var_level]][[metric]] <- model_summary
              
            } # End loop level 4 (diversity/distribution metric) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (LMM model)
            
          } # End loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (LMM model)
          
        } # End loop level 2 (subset variable) << loop level 1 (LMM model)
        
        # # Combine results into a data frame
        # results_df <- bind_rows(da_res_list[[paste0("model_", model_number)]][[method]])
        
      } # End loop level 1 (LMM model)
      
      da_res_df <- bind_rows(rlang::squash(da_res_list)) # `squash` recursively flattens the list
      
    } # End valid formula table check
    
    ### ................................................
    ###
    ### Visualization ----
    ###
    ### ................................................
    if(!is.null(grouping_vars_tcr) & (sum(grouping_vars_tcr == "") < length(grouping_vars_tcr))) {
      
      # Set `subset_vars_tcr` to whatever it was in the config YAML file
      if(subset_vars_tcr_orig_flagged) subset_vars_tcr <- NULL
      # See if we need to do any subsetting prior to graphing
      if(flagVariable(subset_vars_tcr) | "Complete data set" %in% subset_vars_tcr) {
        # If there are no subset variables/one of them is "Complete data set", we will add a column that will act as a dummy subset variable
        # and change `subset_vars_tcr` to be the name of this dummy subset variable
        # This will allow us to use one loop for either case (controls switch 1a or 1b)
        pData(target_data_object)[["Complete data set"]] <- "Complete data set"
        pData(target_data_object)[["Complete data set"]] <- as.factor(pData(target_data_object)[["Complete data set"]])
        subset_vars_tcr <- union(subset_vars_tcr, "Complete data set") # [is.na(subset_vars_tcr)]
        
      } # End control switch 1a (no subset variables) << loop level 1 (model)
      
      # Graph 16S expression levels by group,
      # subsetting if requested
      plot_list <- list()
      
      for(subset_var in subset_vars_tcr) {
        plot_list[[subset_var]] <- list()
        
        # Check if it's NA
        if(is.na(subset_var) || subset_var=="NA") {
          subset_var <- "Complete data set"
        }
        
        # Get the levels
        subset_var_levels <- pData(target_data_object)[[subset_var]] %>% unique
        
        # Loop over the levels and subset by each level
        for(subset_var_level in subset_var_levels) {
          plot_list[[subset_var]][[subset_var_level]] <- list()
          
          pData_sub <- pData(target_data_object) %>% dplyr::filter(!!as.name(subset_var) == subset_var_level)
          for(grouping_var in grouping_vars_tcr) {
            if(grouping_var==subset_var) next
            plot_list[[grouping_var]] <- list()
            
            # Create the `plot_df` data frame for graphing from `pData_sub`
            # We will select only the columns we need
            # and melt into long form so we can facet by metric type
            plot_df <- pData_sub %>% 
              dplyr::select(c(grouping_var, metric_names)) %>%
              reshape2::melt(id.vars = grouping_var,
                             variable.name = "metric")
            
            # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
            excluded_levels <- excluded_levels_list[[grouping_var]]
            if(!flagVariable(excluded_levels)) plot_df <- plot_df %>% dplyr::filter(!(!!as.name(grouping_var) %in% excluded_levels))
            
            # Convert the current grouping variable to factor
            plot_df[[grouping_var]] <- as.factor(plot_df[[grouping_var]])
            
            # Make sure we have enough colors
            n_colors <- plot_df[[grouping_var]] %>% unique %>% length
            mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
            
            # Graph
            plot_title <- ifelse(subset_var == "Complete data set", "", glue::glue("Subset variable: {subset_var} | level: {subset_var_level}"))
            plot <- plot_df %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = value, fill = !!as.name(grouping_var))) + 
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
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) + 
              facet_wrap(c("metric"), ncol = 2, scales = "free") +
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
                   y = paste0(""))
            
            # If `add_boxplot_pvals_tcr` is TRUE, 
            # build p-value data frame and add p-values to plot
            # using `stat_pvalue_manual()`
            if(add_boxplot_pvals_tcr) {
              # Check if the current `grouping_var` is in `da_res_df` (the data frame with all the differential LMM results)
              if(grouping_var %in% da_res_df$fixed_effect) {
                # Create the p-value data frame
                # p-value data frame: group1, group2, p, y.position, p.label 
                pvals_df <- da_res_df %>% 
                  dplyr::filter(fixed_effect == grouping_var & subset_var == !!subset_var & subset_var_level == !!subset_var_level) %>%
                  dplyr::select(baseline, term, p.value, metric) %>% 
                  dplyr::rename(group1 = baseline, group2 = term, p = p.value) %>% 
                  # Because the LMM or whatever reformats values with special characters
                  # by enclosing them in parentheses ("()"), we need to strip those parentheses
                  # out of the values, otherwise it will cause weird stuff to happen with the graphing
                  dplyr::mutate(group1 = group1 %>% regexPipes::gsub("^\\(", "") %>% regexPipes::gsub("\\)$", ""),
                                group2 = group2 %>% regexPipes::gsub("^\\(", "") %>% regexPipes::gsub("\\)$", ""))
                
                # Calculate y.position: slightly above the highest point per facet
                tops <- plot_df %>% # `plot_df` should already have only the samples with the correct `subset_var` and `subset_var_level`
                  group_by(metric) %>% # Group by metric type
                  summarise(y.position = max(value, na.rm = TRUE) * 1.05, .groups = "drop") %>%
                  as.data.frame
                
                # Add `y.position` to `pvals_df`
                # The space between brackets should be ~ 10% the range of the points
                ranges <- plot_df %>% group_by(metric) %>% summarise(range = range(value)) %>% summarise(range = diff(range)) %>% as.data.frame
                bracket_spacing <- 0.15 * ranges[,2]; names(bracket_spacing) <- ranges[,1]
                highest_bracket <- bracket_spacing * ((pvals_df %>% nrow()) / 4 - 1)
                # Calculate the y-positions of the brackets
                y.position <- c()
                for(metric_name in metric_names) {
                  y.position <- c(y.position, (tops %>% dplyr::filter(metric==metric_name) %>% .[1,2]) + seq(from = 0, to = highest_bracket %>% .[names(.)==metric_name], by = bracket_spacing %>% .[names(.)==metric_name]))
                }
                pvals_df[["y.position"]] <- y.position
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
            file_name <- glue::glue("TCR_diversity_plot_{subset_var}_{subset_var_level}_{grouping_var}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 8, height = 6, units = "in")
            }
            
          } # End grouping variables `for()` loop
        } # End subset variable levels `for()` loop
      } # End subset variables `for()` loop
    } # End check for grouping variables
    
    ## ................................................
    ##
    ### Finalization ----
    ##
    ## ................................................
    # Add the main module data object back to the list
    target_data_object_list[[module_tcr]] <- target_data_object
    
    # Now update the pData for the rest of the modules
    for(module_i in names(target_data_object_list)) {
      if(module_i == module_tcr) next
      diversity_cols <- setdiff(colnames(pData(target_data_object)), colnames(pData(target_data_object_list[[module_i]])))
      pData(target_data_object_list[[module_i]]) <- cbind(pData(target_data_object_list[[module_i]]), pData(target_data_object)[,diversity_cols])
    }
    
  } # End if() there are TCR probes

} # End if() the TCR module is provided

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
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_TCR-analysis.rds"))
# Export the raw plots as RDS file
if(exists("plot_list")) plot_list %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_plots-list.rds"))
# Export the LMM results
if(exists("da_res_df")) da_res_df %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_LMM-results.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_tcr_analysis.RData"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)