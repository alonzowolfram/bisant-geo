## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Setting up for differential expression analysis")

## Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
target_data_object_list <- readRDS(cl_args[5])
# Set `main_module` if not set already
modules <- names(target_data_object_list)
if(flagVariable(main_module)) main_module <- modules[1]
rm(modules)
# We'll only need the main module for this one
target_data_object <- target_data_object_list[[main_module]]

# Get the test variables
test_vars <- c()
for(formula in lmm_formulae_de) {
  test_vars <- union(test_vars, extractVariables(formula))
}
# Convert test variables and subset variables to factors
for(test_var in test_vars) {
  pData(target_data_object)[[test_var]] <- pData(target_data_object)[[test_var]] %>% as.factor
}
for(subset_var in subset_vars_de) {
  if(!is.null(subset_var) && !is.na(subset_var) && subset_var != "NA") pData(target_data_object)[[subset_var]] <- pData(target_data_object)[[subset_var]] %>% as.factor
}

# # Create the table of LMM parameter combinations
# param_combos <- cbind(random_slope, test_vars, random_intercept_vars, random_slope_vars) %>% as.data.frame %>% dplyr::mutate(`Model number` = 1:nrow(.)) # , random_slope_vars
# colnames(param_combos)[1:4] <- c("Random slope", "Test variable", "Random intercept variable", "Random slope variable") # , "Random slope variable"

# Get the negative probes. This code was used in qc_probes.R, and
# maybe in the future, we will have a way to remove this redundant code. Maybe by saving the negative probes
# into the NanoStringGeoMxSet object?
negativeProbefData <- subset(fData(target_data_object), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

message("Finished setup for differential expression analysis")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Differential expression analysis via LMM ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Run LMM:
# formula follows conventions defined by the lme4 package
# When running LMM, mixedModelDE seems to have issues with variable names with spaces etc., even if enclosed in backticks ("``")

# Now we can loop over all the param combos and run LMM on each one
message("Performing differential expression analysis")

results2 <- data.frame()
model_number <- 1
for (lmm_formula in lmm_formulae_de) {
  message(glue::glue("Working on model {lmm_formula}"))
  
  # Create the model formula
  model_formula <- as.formula(lmm_formula)
  # Extract the grouping variable (by convention, it will be the first fixed effect in the model)
  groupVar <- extractFirstFixedEffect(model_formula)
  
  # Handle missing subset variables
  if (all(is.na(subset_vars_de)) || all(subset_vars_de == "NA", na.rm = TRUE)) {
    pData(target_data_object)$`All observations` <- factor("DummyLevel")
    subset_vars_de <- "All observations"
  }
  
  for (subset_var in subset_vars_de) {
    message(glue::glue("Working on subset variable {subset_var} for model {lmm_formula}"))
    
    if (subset_var %in% c("NA", NA)) {
      pData(target_data_object)$`All observations` <- factor("DummyLevel")
      subset_var <- "All observations"
    }
    
    model_vars <- extractModelComponents(model_formula) %>% unlist
    if (subset_var %in% model_vars) {
      message(glue::glue("Skipping subset variable {subset_var} due to overlap with one or more fixed-effect or random variables in the model formula"))
      next
    }
    
    subset_levels <- levels(factor(pData(target_data_object)[[subset_var]]))
    # If subset_var_de_levels_manual is set, filter subset_var_levels to include only those values
    subset_var_de_levels_manual_i <- subset_var_de_levels_manual[[subset_var]]
    if(sum(is.na(subset_var_de_levels_manual_i)) < length(subset_var_de_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
      if(sum(subset_var_de_levels_manual_i %in% subset_levels) < 1) {
        warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of {subset_var} will be used instead"))
      } else { # At least one subset_var_level_manual_i is an actual level of the current subset variable
        subset_levels <- subset_levels %>% .[. %in% subset_var_de_levels_manual_i]
      }
    }
    
    for (subset_level in subset_levels) {
      message(glue::glue("Working on level {subset_level} of {subset_var} for model {lmm_formula}"))
      
      ind <- pData(target_data_object)[[subset_var]] == subset_level
      
      for (norm_method in normalization_methods) {
        message(glue::glue("Working on normalization {norm_method} for level {subset_level} of {subset_var} for model {lmm_formula}"))
        skip_to_next <- FALSE
        
        set.seed(random_seed)
        CV_dat <- assayDataApply(target_data_object, elt = paste0("log_", norm_method), MARGIN = 1, calc_CV) %>% na.omit()
        
        if (is.null(cv_cutoff) || !is.finite(cv_cutoff) || !(0 <= cv_cutoff & cv_cutoff <= 1)) {
          top_cv_genes <- names(CV_dat)
        } else {
          top_cv_genes <- names(CV_dat[CV_dat > quantile(CV_dat, probs = cv_cutoff)])
        }
        
        mean_exprs_cutoff_genes <- rowMeans(target_data_object@assayData[[paste0("log_", norm_method)]]) > 1
        genes_pass_cutoffs <- intersect(top_cv_genes, names(mean_exprs_cutoff_genes[mean_exprs_cutoff_genes]))
        
        tdo_sub <- subset(target_data_object, TargetName %in% genes_pass_cutoffs, ind)
        
        # Fit the model.
        mixedOutmc <- tryCatch(
          mixedModelDE(tdo_sub, 
                       elt = paste0("log_", norm_method), 
                       modelFormula = model_formula, 
                       groupVar = groupVar, 
                       nCores = parallel::detectCores(), multiCore = TRUE),
          error = function(e) {
            warning(glue::glue("Error in model {lmm_formula} - {norm_method} - {subset_level} - {subset_var}: {e$message}"))
            skip_to_next <<- TRUE
          }
        )
        
        if (skip_to_next) next
        
        r_test <- do.call(rbind, mixedOutmc["lsmeans", ]) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column(var = "Contrast") %>% 
          mutate(
            Contrast = Contrast %>% pipe.gsub("\\.[[:digit:]]+$", "") %>% pipe.gsub("\\.{3}", " \\- "),
            Gene = rep(colnames(mixedOutmc), each = n()/ncol(mixedOutmc)),
            `Subset variable` = subset_var,
            `Subset level` = subset_level,
            `Normalization method` = normalization_names[names(normalization_names) == norm_method],
            FDR = p.adjust(`Pr(>|t|)`, method = "fdr"),
            `Contrast variable` = groupVar,
            `Model` = as.character(lmm_formula)
            # ,`Model number` = i
        ) %>% 
          relocate(`Contrast variable`, .before = Contrast) %>%
          relocate(Gene, .before = Estimate)
        
        results2 <- bind_rows(results2, r_test)
      }
    }
  }
  
  if ("All observations" %in% colnames(pData(target_data_object))) {
    pData(target_data_object) <- dplyr::select(pData(target_data_object), -`All observations`)
  }
  
  # colnames(pData(target_data_object)) <- orig_var_names
  if (any(subset_vars_de == "All observations")) subset_vars_de <- NA
}
message("Differential expression analysis completed")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Volcano plots ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Graphing differentially expressed genes")

## ................................................
##
### Selection of genes for labels ----
##
## ................................................
# Order genes for convenience:
results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
top_g_list <- list()

models <- results2$`Model` %>% unique
model_number <- 1
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it "subset_vars_de" because we already have a variable by that name ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  # If subset_var_de_levels_manual is set, filter subset_var_levels to include only those values
  subset_var_de_levels_manual_i <- subset_var_de_levels_manual[[subset_var]]
  if(sum(is.na(subset_var_de_levels_manual_i)) < length(subset_var_de_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
    if(sum(subset_var_de_levels_manual_i %in% subset_var_levels) < 1) {
      warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of {subset_var} will be used instead"))
    } else { # At least one subset_var_level_manual_i is an actual level of the current subset variable
      subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_de_levels_manual_i]
    }
  }
  
  for(subset_var_level in subset_var_levels) {
    
    for(model in models) {
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & 
                                                `Subset level`==subset_var_level & 
                                                `Model`==as.character(model)) %>% .$Contrast %>% unique
      for(contrast in contrasts) {
        
        for(normalization_method in unique(results2$`Normalization method`)) {
          top_g <- c()
          ind <- results2$`Subset variable`==subset_var & 
            results2$`Subset level`==subset_var_level & 
            results2$`Model`==as.character(model) & 
            results2$Contrast==contrast & 
            results2$`Normalization method`==normalization_method
          
          # Populate `top_g`
          top_g <- c(top_g,
                     results2[ind, 'Gene'][
                       order(results2[ind, 'invert_P'], decreasing = TRUE)[1:n_top_genes]],
                     results2[ind, 'Gene'][
                       order(results2[ind, 'invert_P'], decreasing = FALSE)[1:n_top_genes]]) %>% unique
          top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]] <- top_g
          
        }
      }
      
      # Increment `model_number`
      model_number <- model_number + 1
      
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
# Set the significance colors
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
# Initialize the list to hold the plots
plot_list_diff_exprs <- list()
model_list <- list() # Hold all the model formulae so we can access them when arranging the plots into grids
models <- results2$`Model` %>% unique
model_number <- 1
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it subset_vars_de because we already have a variable by that name ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model in models) {
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model`==model) %>% .$Contrast %>% unique
      for(normalization_method in unique(results2$`Normalization method`)) {
        
        for(contrast in contrasts) {
          # Get the parameters for this experiment
          test_var <- results2 %>% dplyr::filter(`Model`==model) %>% dplyr::select(`Contrast variable`) %>% unlist %>% .[1]
          test_var_lv_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1] #pData(target_data_object)[[test_var]] %>% levels %>% .[1]
          test_var_lv_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2] #pData(target_data_object)[[test_var]] %>% levels %>% .[2]
          if(subset_var=="NA" || is.na(subset_var) || subset_var=="All observations") {
            subset_by <- ""
          } else {
            subset_by <- glue::glue("| Subset variable: {subset_var}, level: {subset_var_level}")
          }
          
          # Get the top DE genes for this experiment
          # top_genes <- c()
          # for(subset_var_level in subset_var_levels) {
          #   top_genes <- c(top_genes, top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]])
          # }
          # top_genes <- unique(top_genes)
          # Above code is a remnant from when we faceted using facet_by() and had to include all top genes
          # Now that we manually create a grid, we can just include the appropriate top genes
          top_genes <- top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]][[normalization_method]]
          
          # Set up the data to plot
          n <- length(contrasts)
          nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # Used for the font size of the volcano plot labels
          
          results2_sub <- results2 %>% dplyr::filter(
            `Subset variable`==subset_var &
              `Subset level`==subset_var_level &
              `Model`==model &
              Contrast==contrast & 
              `Normalization method`==normalization_method
          )
          # Make the ggplot
          plot <- ggplot(results2_sub,
                         aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                             color = Color, label = Gene)) +
            geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
            geom_hline(yintercept = -log10(0.05), lty = "dashed") +
            geom_point() +
            scale_color_manual(values = signif_cols) +
            labs(x = glue::glue("Enriched in {test_var_lv_2} ←  log2(FC) → Enriched in {test_var_lv_1}"),
                 y = "", # Significance, -log10(P)
                 # title = paste0("DE genes", subset_by, " \nTest variable: ", test_var, random_slope_status),
                 color = "Significance") +
            scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
            geom_text_repel(data = subset(results2_sub, Gene %in% top_genes), # & FDR < 0.001 # The way we have the graphing for the DEGs set up, it will label all genes for a given subset variable, across all values of that variable. This is because we use the values of the subset variable to facet, and AFAIK, ggplot2 doesn't have a way to exclude labels by the variable that's being faceted by. This may change
                            size = (4 / nCol), point.padding = 0.15, color = "black",
                            min.segment.length = .1, box.padding = .2, lwd = 2,
                            max.overlaps = 50) +
            theme_bw(base_size = 16) +
            theme(legend.position = "bottom",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()) # "bottom"
          
          # Add to the list
          plot_list_diff_exprs[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[normalization_method]][[contrast]] <- plot
          
        }
      }
      
      model_list[[paste0("model_", model_number)]] <- as.character(model)
      model_number <- model_number + 1
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
          # Stop right there
          warning(glue::glue("The combination of {subset_var} - {subset_var_level} - {model_num} -  {normalization_method} has {n} volcano plots, too many for graphing. Please graph these manually. Skipping to the next list of graphs"))
          rm(p_list)
          gc()
          next
        }
        nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1
        
        # Set the scaling factors for label and legend size
        sqrt_n_col <- sqrt(nCol)
        scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 2), 1) # Number of rows in current grid / 2 (base number)
        
        # Create new dummy plot, scale legend accordingly, and then extract legend
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
        model <- model_list[[model_num]] %>% as.character
        test_var <- results2 %>% dplyr::filter(Model == model) %>% .$`Contrast variable` %>% unique %>% .[1]
        if(subset_var=="NA" || is.na(subset_var) || subset_var=="All observations") {
          subset_by <- ""
        } else {
          subset_by <-glue::glue("| Subset variable: {subset_var}, level: {subset_var_level}")
        }
        
        # Arrange plots in p_list onto a grid
        plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
        plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Significance, -log10(P)", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = scaling_factor)),
                                                           bottom = grid::textGrob("", gp = grid::gpar(cex = scaling_factor)),
                                                           top = grid::textGrob(glue::glue("DE genes {subset_by} \nTest (contrast) variable: {test_var} | {model} \nNormalization method: {normalization_method}"), 
                                                                                gp = grid::gpar(cex = scaling_factor)))
        
        # Add back in the legend we extracted earlier
        plot_grid2 <- grid.arrange(plot_grid, legend, ncol = 1, heights=c(10, 1))
        
        # Add to the plot list
        plot_grid2 <- ggplotify::as.ggplot(plot_grid2)
        plot_list_diff_exprs_grid[[subset_var]][[subset_var_level]][[model_num]][[normalization_method]] <- plot_grid2
        # Save plot to disk
        plot_grid2
        file_name <- glue::glue("DE_volcano_plot_grid_{subset_var}_{subset_var_level}_{model_num}_{normalization_method}")
        for(file_type in output_plot_file_types) {
          ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
        }
        
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
message("Exporting differential expression analysis results")

# Export tables of DE genes to CSV.
results2 %>% write.csv(paste0(output_dir_tabular, "LMM-differential-expression_results.csv")) 
# Export graphs.
saveRDS(plot_list_diff_exprs, paste0(output_dir_rdata, "LMM-DEG_volcano-plots.rds"))
saveRDS(plot_list_diff_exprs_grid, paste0(output_dir_rdata, "LMM-DEG_volcano-plot_grids.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_differential_expression_analysis.RData"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)

message("Differential expression analysis results successfully exported")