## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Setting up for differential expression analysis.")

## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
target_data_object_list <- readRDS(cl_args[5])
# We'll only need the main module for this one.
target_data_object <- target_data_object_list[[main_module]]

# Convert test variables and subset variables to factors.
for(test_var in test_vars) {
  pData(target_data_object)[[test_var]] <- pData(target_data_object)[[test_var]] %>% as.factor
}
for(subset_var in subset_vars) {
  if(!is.null(subset_var) && !is.na(subset_var) && subset_var != "NA") pData(target_data_object)[[subset_var]] <- pData(target_data_object)[[subset_var]] %>% as.factor
}

# Create the table of LMM parameter combinations.
param_combos <- cbind(random_slope, test_vars, random_intercept_vars) %>% as.data.frame %>% dplyr::mutate(`Model number` = 1:nrow(.)) # , random_slope_vars
colnames(param_combos)[1:3] <- c("Random slope", "Test variable", "Random intercept variable") # , "Random slope variable"

# Get the negative probes. This code was used in qc_probes.R, and
# maybe in the future, we will have a way to remove this redundant code. Maybe by saving the negative probes
# into the NanoStringGeoMxSet object?
negativeProbefData <- subset(fData(target_data_object), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

message("Finished setup for differential expression analysis.")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Differential expression analysis via LMM ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Run LMM:
# formula follows conventions defined by the lme4 package.
# When running LMM, mixedModelDE seems to have issues with variable names with spaces etc., even if enclosed in backticks (``). 
# So we're going to rename the variables of interest to something without spaces and other special characters.
# But first, we'll have to save the original variable names so every time we rename the variables in the for() loop,.
# we can reset them afterwards.
orig_var_names <- colnames(pData(target_data_object))
results2 <- c()

# Loop levels (outermost to innermost):
# 1) Model (param_combos)
# 2) Subset variable (subset_vars)

# Now we can loop over all the param combos and run LMM on each one.
# Loop level 1: model.
message("Performing differential expression analysis.")
for(i in 1:nrow(param_combos)) {
  message(paste0("Working on model #", i, "."))
  
  # Get the params for this experiment. 
  random_slope_i <- param_combos[i,1]
  test_var <- param_combos[i,2]
  random_intercept_var <- param_combos[i,3]
  # random_slope_vars <- param_combos[i,4]
  
  # Reset the variable names.
  colnames(pData(target_data_object)) <- orig_var_names
  
  # Rename the test variable (test_var) and the random intercept variable (random_intercept_var).
  pData(target_data_object) <- pData(target_data_object) %>%
    dplyr::rename(TestVar = !!as.name(test_var), RandomInterceptVar = !!as.name(random_intercept_var))
  
  # control switch 1 (subset variables exist or not) << loop level 1 (model).
  if(sum(is.na(subset_vars)) == length(subset_vars) || sum(subset_vars=="NA", na.rm = T) == length(subset_vars[!is.na(subset_vars)])) {
    # control switch 1a: no subset variables << loop level 1 (model).
    
    # Since there are no subset variables, we will add a column that will act as a dummy subset variable
    # and change subset_vars to be the name of this dummy subset variable.
    # This will allow us to use one loop for either case (controls switch 1a or 1b).
    pData(target_data_object)[["All observations"]] <- "DummyLevel"
    pData(target_data_object)[["All observations"]] <- as.factor(pData(target_data_object)[["All observations"]])
    subset_vars <- c("All observations")
    
  } # End control switch 1a (no subset variables) << loop level 1 (model).
  
  # loop level 2 (subset variable) << loop level 1 (model)
  for(subset_var in subset_vars) {
    message(paste0("Working on subset variable ", subset_var, " for model #", i, "."))
    
    # Check if the current subset_var is NA.
    if(subset_var == "NA" || is.na(subset_var)) {
      # Add a column that will act as a dummy subset variable
      # and change subset_vars to be the name of this dummy subset variable.
      # This will allow us to use one loop for either case.
      pData(target_data_object)[["All observations"]] <- "DummyLevel"
      pData(target_data_object)[["All observations"]] <- as.factor(pData(target_data_object)[["All observations"]])
      subset_var <- "All observations"
    }
    
    # Check that the current subset_var is not the same as either test_var or random_intercept_var.
    if(test_var==subset_var || random_intercept_var==subset_var) {
      message(paste0("The current subset variable, ", subset_var, " is the same as either your test (contrast) variable or random intercept variable. Skipping this subset variable - contrast/random intercept variable combination."))
      next
    }
      
    # Get the levels of the current subset_var.
    subset_var_levels <- pData(target_data_object)[[subset_var]] %>% as.factor %>% levels # as.factor needed because it might be a character vector.
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      message(paste0("Working on level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
      
      # Get all the samples belonging to the current subset_var_level.
      ind <- pData(target_data_object)[[subset_var]] == subset_var_level
      
      # Set the model formula based on whether we have a random slope or not. 
      model_formula <- ifelse(random_slope_i %in% c("no", "FALSE"), "~ TestVar + (1 | RandomInterceptVar)", "~ TestVar + (1 + TestVar | RandomInterceptVar)") %>% as.formula
        
      # loop level 4 (current normalization method) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
      for(norm_method in normalization_methods) {
        message(paste0("Working on normalization method ", norm_method, " for level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
        # Error catching: https://stackoverflow.com/a/55937737/23532435
        skip_to_next <- FALSE
        
        # # Set up the data. 
        # # Expresssion: rows are GENES and columns are SAMPLES. 
        # # Also, remove negative probes.
        # exprs <- target_data_object@assayData[[paste0("log_", norm_method)]] %>% .[!(rownames(.) %in% neg_probes),]
        # metadata <- sData(target_data_object)
        
        # Set the random seed.
        set.seed(random_seed)
        # Calculate coefficient of variance for each gene.
        CV_dat <- assayDataApply(target_data_object,
                                 elt = paste0("log_", norm_method), MARGIN = 1, calc_CV) %>% .[!is.na(.)]
        # Keep only the genes with the highest CVs
        # 2024/12/16: and genes with mean expression (across samples) > 1.
        # Genes that meet CV cutoff.
        mean_cv <- mean(CV_dat, na.rm = T)
        sd_cv <- sd(CV_dat, na.rm = T)
        if(is.null(cv_cutoff) || cv_cutoff == "" || !is.finite(cv_cutoff) || !(0 <= cv_cutoff & cv_cutoff <= 1)) {
          message("cv_cutoff for genes either has not been provided or is not a numeric value in [0,1]. Using all genes for differential expression.")
          top_cv_genes <- CV_dat %>% names
        } else {
          message(paste0("Using cv_cutoff of ", cv_cutoff, " for differential expression."))
          nth_percentile_cv <- quantile(CV_dat, probs = cv_cutoff)
          top_cv_genes <- CV_dat %>% .[. > nth_percentile_cv] %>% names
          # top_cv_genes <- CV_dat %>% .[. > (mean_cv + cv_cutoff*sd_cv)] %>% names
        }
        # Genes that meet mean expression cutoff.
        mean_exprs_cutoff <- 1
        mean_exprs_cutoff_genes <- which(rowMeans(target_data_object@assayData[[paste0("log_", norm_method)]]) > mean_exprs_cutoff) %>% names
        # Union of the two. 
        # EDIT 2025/01/06: Should be intersection of the two.
        genes_pass_cutoffs <- intersect(top_cv_genes, mean_exprs_cutoff_genes) # union(top_cv_genes, mean_exprs_cutoff_genes)
        
        # Generate the model.
        tdo_sub <- subset(target_data_object, TargetName %in% genes_pass_cutoffs, ind)
        mixedOutmc <- tryCatch(mixedModelDE(tdo_sub, # Error handling needed because sometimes there might not be enough samples to perform DE.
                                   elt = paste0("log_", norm_method),
                                   modelFormula = model_formula,
                                   groupVar = "TestVar",
                                   nCores = parallel::detectCores(),
                                   multiCore = TRUE),
                               error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next) {
          warning(paste0("An error occurred when trying to perform differential expression for normalization method ", norm_method, " for level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
          next
        }
        
        # Format results as data.frame.
        r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
        tests <- rownames(r_test)
        r_test <- as.data.frame(r_test)
        r_test$Contrast <- tests
        
        # Use lapply in case you have multiple levels of your test factor to
        # correctly associate gene name with its row in the results table.
        r_test$Gene <- 
          unlist(lapply(colnames(mixedOutmc),
                        rep, nrow(mixedOutmc["lsmeans", ][[1]])))
        r_test$`Subset variable` <- subset_var
        r_test$`Subset level` <- subset_var_level
        r_test$`Normalization method` <- normalization_names[names(normalization_names)==norm_method]
        r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
        r_test <- r_test[, c("Gene", "Subset variable", "Subset level", "Normalization method", "Contrast", "Estimate", 
                             "Pr(>|t|)", "FDR")]
        r_test$`Contrast variable` <- test_var
        r_test$`Model number` <- i
        r_test <- r_test %>% dplyr::relocate(`Contrast variable`, .before = Contrast)
        results2 <- rbind(results2, r_test)
        
      } # End loop level 4: normalization method. 
      
    } # End loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model).
    
    # Now check if we have any manual subset levels.
    subset_var_levels_manual_i <- subset_var_levels_manual[[subset_var]]
    if(!is.null(subset_var_levels_manual_i)) {
      if(sum(is.na(subset_var_levels_manual_i)) < length(subset_var_levels_manual_i)) {
        message(paste0("Working on levels ", paste(subset_var_levels_manual_i, collapse = ", "), " for subset variable ", subset_var, " for model #", i, "."))
        
        # Get all the samples belonging to the current subset_var_level.
        ind <- pData(target_data_object)[[subset_var]] %in% subset_var_levels_manual_i
        
        # Set the model formula based on whether we have a random slope or not. 
        model_formula <- ifelse(random_slope_i %in% c("no", "FALSE"), "~ TestVar + (1 | RandomInterceptVar)", "~ TestVar + (1 + TestVar | RandomInterceptVar)") %>% as.formula
        
        # loop level 4 (current normalization method) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
        for(norm_method in normalization_methods) {
          message(paste0("Working on normalization method ", norm_method, " for levels ", paste(subset_var_levels_manual_i, collapse = ", "), " for subset variable ", subset_var, " for model #", i, "."))
          # Error catching: https://stackoverflow.com/a/55937737/23532435
          skip_to_next <- FALSE
          
          # # Set up the data. 
          # # Expresssion: rows are GENES and columns are SAMPLES. 
          # # Also, remove negative probes.
          # exprs <- target_data_object@assayData[[paste0("log_", norm_method)]] %>% .[!(rownames(.) %in% neg_probes),]
          # metadata <- sData(target_data_object)
          
          # Set the random seed.
          set.seed(random_seed)
          # Calculate coefficient of variance for each gene.
          CV_dat <- assayDataApply(target_data_object,
                                   elt = paste0("log_", norm_method), MARGIN = 1, calc_CV) %>% .[!is.na(.)]
          # Keep only the genes with the highest CVs
          # 2024/12/16: and genes with mean expression (across samples) > 1.
          # Genes that meet CV cutoff.
          mean_cv <- mean(CV_dat, na.rm = T)
          sd_cv <- sd(CV_dat, na.rm = T)
          if(is.null(cv_cutoff) || cv_cutoff == "" || !is.finite(cv_cutoff) || !(0 <= cv_cutoff & cv_cutoff <= 1)) {
            message("cv_cutoff for genes either has not been provided or is not a numeric value in [0,1]. Using all genes for differential expression.")
            top_cv_genes <- CV_dat %>% names
          } else {
            message(paste0("Using cv_cutoff of ", cv_cutoff, " for differential expression."))
            nth_percentile_cv <- quantile(CV_dat, probs = cv_cutoff)
            top_cv_genes <- CV_dat %>% .[. > nth_percentile_cv] %>% names
            # top_cv_genes <- CV_dat %>% .[. > (mean_cv + cv_cutoff*sd_cv)] %>% names
          }
          # Genes that meet mean expression cutoff.
          mean_exprs_cutoff <- 1
          mean_exprs_cutoff_genes <- which(rowMeans(target_data_object@assayData[[paste0("log_", norm_method)]]) > mean_exprs_cutoff) %>% names
          # Union of the two. 
          # EDIT 2025/01/06: Should be intersection of the two.
          genes_pass_cutoffs <- intersect(top_cv_genes, mean_exprs_cutoff_genes) # union(top_cv_genes, mean_exprs_cutoff_genes)
          
          # Generate the model.
          tdo_sub <- subset(target_data_object, TargetName %in% genes_pass_cutoffs, ind)
          mixedOutmc <- tryCatch(mixedModelDE(tdo_sub, # Error handling needed because sometimes there might not be enough samples to perform DE.
                                              elt = paste0("log_", norm_method),
                                              modelFormula = model_formula,
                                              groupVar = "TestVar",
                                              nCores = parallel::detectCores(),
                                              multiCore = TRUE),
                                 error = function(e) {skip_to_next <<- TRUE})
          if(skip_to_next) {
            warning(paste0("An error occurred when trying to perform differential expression for normalization method ", norm_method, " for level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
            next
          }
          
          # Format results as data.frame.
          r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
          tests <- rownames(r_test)
          r_test <- as.data.frame(r_test)
          r_test$Contrast <- tests
          
          # Use lapply in case you have multiple levels of your test factor to
          # correctly associate gene name with its row in the results table.
          r_test$Gene <- 
            unlist(lapply(colnames(mixedOutmc),
                          rep, nrow(mixedOutmc["lsmeans", ][[1]])))
          r_test$`Subset variable` <- subset_var
          r_test$`Subset level` <- paste(subset_var_levels_manual_i, collapse = ",")
          r_test$`Normalization method` <- normalization_names[names(normalization_names)==norm_method]
          r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
          r_test <- r_test[, c("Gene", "Subset variable", "Subset level", "Normalization method", "Contrast", "Estimate", 
                               "Pr(>|t|)", "FDR")]
          r_test$`Contrast variable` <- test_var
          r_test$`Model number` <- i
          r_test <- r_test %>% dplyr::relocate(`Contrast variable`, .before = Contrast)
          results2 <- rbind(results2, r_test)
          
        } # End loop level 4: normalization method. 
      }
    } # End if there are manual subset levels.
    
  } # End loop level 2 (subset variable) << loop level 1 (model).
  
  # Reset the variable names. 
  # First check if there's an additional column added if there are no subset variables (see control switch 1a.)
  # If there is, remove it. 
  if("All observations" %in% colnames(pData(target_data_object))) pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-`All observations`)
  # Then reset the variable names. 
  colnames(pData(target_data_object)) <- orig_var_names 
  
  # If we changed subset_vars to "All observations", change it back to NA.
  if(sum(subset_vars=="All observations", na.rm = T) == 1) subset_vars <- NA
  
} # End loop level 1: model.
message("Differential expression analysis completed.")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Volcano plots ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Graphing differentially expressed genes.")

## ................................................
##
### Selection of genes for labels ----
##
## ................................................
# Order genes for convenience:
results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
top_g_list <- list()

model_numbers <- results2$`Model number` %>% unique
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it subset_vars because we already have a variable by that name. ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model_number in model_numbers) {
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & 
                                                `Subset level`==subset_var_level & 
                                                `Model number`==model_number) %>% .$Contrast %>% unique
      for(contrast in contrasts) {
        
        for(normalization_method in unique(results2$`Normalization method`)) {
          top_g <- c()
          ind <- results2$`Subset variable`==subset_var & 
            results2$`Subset level`==subset_var_level & 
            results2$`Model number`==model_number & 
            results2$Contrast==contrast & 
            results2$`Normalization method`==normalization_method
          
          # Populate top_g.
          top_g <- c(top_g,
                     results2[ind, 'Gene'][
                       order(results2[ind, 'invert_P'], decreasing = TRUE)[1:n_top_genes]],
                     results2[ind, 'Gene'][
                       order(results2[ind, 'invert_P'], decreasing = FALSE)[1:n_top_genes]]) %>% unique
          top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]] <- top_g
          
        }
      }
    }
  }
}
results2 <- results2[, -1*ncol(results2)] # remove invert_P from matrix

## ................................................
##
### Graphing the volcano plots ----
##
## ................................................
# Categorize results2 based on P-value & FDR for plotting
results2$Color <- "NS or FC < 0.5"
results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results2$Color[results2$FDR < 0.25] <- "FDR < 0.25"
results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
results2$Color <- factor(results2$Color,
                         levels = c("NS or FC < 0.5", 
                                    "P < 0.05",
                                    "FDR < 0.25",
                                    "FDR < 0.05", 
                                    "FDR < 0.001"))
# Set the significance colors.
signif_cols <- c("dodgerblue",
                  "lightblue",
                  "orange2",
                  "khaki1",
                  "grey")
names(signif_cols) <- c("FDR < 0.001", 
                        "FDR < 0.05",  
                        "FDR < 0.25",
                        "P < 0.05", 
                        "NS or FC < 0.5")

# Graph results2
# Initialize the list to hold the plots.
plot_list_diff_exprs <- list()
model_numbers <- results2$`Model number` %>% unique
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it subset_vars because we already have a variable by that name. ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model_number in model_numbers) {
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model number`==model_number) %>% .$Contrast %>% unique
      for(normalization_method in unique(results2$`Normalization method`)) {
        
        for(contrast in contrasts) {
          # Get the parameters for this experiment.
          test_var <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Test variable`) %>% unlist %>% .[1]
          random_slope_i <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Random slope`) %>% unlist %>% .[1]
          random_slope_status <- ifelse(random_slope_i %in% c("no", "FALSE"), " | No random slope", " | With random slope")
          test_var_lv_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1] #pData(target_data_object)[[test_var]] %>% levels %>% .[1]
          test_var_lv_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2] #pData(target_data_object)[[test_var]] %>% levels %>% .[2]
          if(subset_var=="NA" || is.na(subset_var) || subset_var=="All observations") {
            subset_by <- ""
          } else {
            subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
          }
          
          # Get the top DE genes for this experiment.
          # top_genes <- c()
          # for(subset_var_level in subset_var_levels) {
          #   top_genes <- c(top_genes, top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]])
          # }
          # top_genes <- unique(top_genes)
          # Above code is a remnant from when we faceted using facet_by() and had to include all top genes.
          # Now that we manually create a grid, we can just include the appropriate top genes.
          top_genes <- top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]]
          
          # Plot.
          n <- length(contrasts)
          nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # Used for the font size of the volcano plot labels.
          
          results2_sub <- results2 %>% dplyr::filter(
            `Subset variable`==subset_var &
              `Subset level`==subset_var_level &
              `Model number`==model_number &
              Contrast==contrast & 
              `Normalization method`==normalization_method
          )
          plot <- ggplot(results2_sub,
                         aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                             color = Color, label = Gene)) +
            geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
            geom_hline(yintercept = -log10(0.05), lty = "dashed") +
            geom_point() +
            scale_color_manual(values = signif_cols) +
            labs(x = paste0("Enriched in ", test_var_lv_2, " <- log2(FC) -> Enriched in ", test_var_lv_1),
                 y = "", # Significance, -log10(P)
                 # title = paste0("DE genes", subset_by, " \nTest variable: ", test_var, random_slope_status),
                 color = "Significance") +
            scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
            geom_text_repel(data = subset(results2_sub, Gene %in% top_genes), # & FDR < 0.001 # The way we have the graphing for the DEGs set up, it will label all genes for a given subset variable, across all values of that variable. This is because we use the values of the subset variable to facet, and AFAIK, ggplot2 doesn't have a way to exclude labels by the variable that's being faceted by. This may change.
                            size = (4 / nCol), point.padding = 0.15, color = "black",
                            min.segment.length = .1, box.padding = .2, lwd = 2,
                            max.overlaps = 50) +
            theme_bw(base_size = 16) +
            theme(legend.position = "bottom",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()) # "bottom"
          
          # Add to the list.
          plot_list_diff_exprs[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[normalization_method]][[contrast]] <- plot
          
        }
      }
    }
  }
}

## ................................................
##
### Arranging volcano plots into grids ----
##
## ................................................
plot_list_diff_exprs_grid <- list()
for(subset_var in names(plot_list_diff_exprs)) {
  for(subset_var_level in names(plot_list_diff_exprs[[subset_var]])) {
    for(model_num in names(plot_list_diff_exprs[[subset_var]][[subset_var_level]])) {
      for(normalization_method in names(plot_list_diff_exprs[[subset_var]][[subset_var_level]][[model_num]])) {
        # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
        # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
        # https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0
        
        p_list <- plot_list_diff_exprs[[subset_var]][[subset_var_level]][[model_num]][[normalization_method]]
        
        n <- length(p_list)
        if(n > 16) {
          # Stop right there.
          warning(paste0("The combination of ", subset_var, " - ", subset_var_level, " - ", model_num, " - ", normalization_method, " has ", n, " volcano plots, too many for graphing. Please graph these manually. Skipping to the next list of graphs."))
          rm(p_list)
          gc()
          next
        }
        nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
        
        # Set the scaling factors for label and legend size.
        sqrt_n_col <- sqrt(nCol)
        scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 2), 1) # Number of rows in current grid / 2 (base number)
        
        # Create new dummy plot, scale legend accordingly, and then extract legend.
        plot_dummy <- p_list[[1]] + 
          scale_color_manual(values = signif_cols,
                             guide = guide_legend(override.aes = list(size = 4 * scaling_factor))
          ) +
          theme(
            # legend.box.background = element_rect(color = "black"),
            legend.title = element_text(size = 14 * scaling_factor),
            legend.key.size = unit(30 * scaling_factor, "points"),
            legend.text = element_text(size = 12 * scaling_factor),
            # legend.key = element_rect(colour = "black"),
            # legend.box.margin = margin(20, 20, 20, 20),
            legend.position = "bottom"
          ) 
        legend <- cowplot::get_plot_component(plot_dummy, 'guide-box-bottom', return_all = TRUE)
        
        # Strip legends from p_list.
        for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}
        
        # Get the model information (test/contrast variable, random slope status, etc.)
        i <- model_num %>% regexPipes::gsub("model_", "")
        random_slope_i <- param_combos[i,1]
        test_var <- param_combos[i,2]
        random_intercept_var <- param_combos[i,3]
        random_slope_status <- ifelse(random_slope_i %in% c("no", "FALSE"), " | No random slope", " | With random slope")
        if(subset_var=="NA" || is.na(subset_var) || subset_var=="All observations") {
          subset_by <- ""
        } else {
          subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
        }
        
        # Arrange plots in p_list onto a grid.
        plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
        plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Significance, -log10(P)", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = scaling_factor)),
                                                           bottom = grid::textGrob("", gp = grid::gpar(cex = scaling_factor)),
                                                           top = grid::textGrob(paste0("DE genes ", subset_by, 
                                                                                       " \nTest (contrast) variable: ", test_var, random_slope_status,
                                                                                       " \nNormalization method: ", normalization_method), 
                                                                                gp = grid::gpar(cex = scaling_factor)))
        
        # Add back in the legend we extracted earlier. 
        plot_grid2 <- grid.arrange(plot_grid, legend, ncol = 1, heights=c(10, 1))
        
        # Save to list.
        plot_list_diff_exprs_grid[[subset_var]][[subset_var_level]][[model_num]][[normalization_method]] <- ggplotify::as.ggplot(plot_grid2)
        
        rm(plot_grid2, p_list)
        gc()
      }
    }
  }
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting differential expression analysis results.")

# Export NanoStringGeoMxSet.
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_differential-expression.rds"))
# Export graphs.
saveRDS(plot_list_diff_exprs, paste0(output_dir_rdata, "LMM-DEG_volcano-plots.rds"))
saveRDS(plot_list_diff_exprs_grid, paste0(output_dir_rdata, "LMM-DEG_volcano-plot_grids.rds"))
# Export tables of DE genes to CSV.
results2 %>% write.csv(paste0(output_dir_tabular, "LMM-differential-expression_results.csv")) 

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)

message("Differential expression analysis results successfully exported.")