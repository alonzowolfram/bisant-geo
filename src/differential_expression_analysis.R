## Source the setup.R file.
source("src/setup.R")

## ---------------------------
# Setup

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

# Create a log2 transform of the data for analysis.
assayDataElement(object = target_data_object, elt = "log_norm", validate = FALSE) <- # Have to set validate to FALSE; otherwise it thinks the dimensions aren't the same. ... 
  assayDataApply(target_data_object, 2, FUN = log, base = 2, elt = normalization_method)

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
# Add to the PowerPoint. 
# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Differential expression analysis (linear mixed models)"), 
                   location = ph_location_label(ph_label = "Title 1"))
# Add the param combos table.
pptx <- pptx %>%
  officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Models created in this run"),
                   location = ph_location_label(ph_label = "Title 1")) %>% 
  officer::ph_with(value = param_combos,
                   location = ph_location_label(ph_label = "Content Placeholder 2"))

## ---------------------------
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
for(i in 1:nrow(param_combos)) {
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
  if(sum(is.na(subset_vars)) == length(subset_vars) || sum(subset_vars=="NA") == length(subset_vars)) {
    # control switch 1a: no subset variables << loop level 1 (model).
    
    # Since there are no subset variables, we will add a column that will act as a dummy subset variable
    # and change subset_vars to be the name of this dummy subset variable.
    # This will allow us to use one loop for either case (controls switch 1a or 1b).
    pData(target_data_object)[["DummySubsetVar"]] <- "DummyLevel"
    pData(target_data_object)[["DummySubsetVar"]] <- as.factor(pData(target_data_object)[["DummySubsetVar"]])
    subset_vars <- c("DummySubsetVar")
    
  } # End control switch 1a (no subset variables) << loop level 1 (model).
  
  # loop level 2 (subset variable) << loop level 1 (model)
  for(subset_var in subset_vars) {
    # Check if the current subset_var is NA.
    if(subset_var == "NA" || is.na(subset_var)) {
      # Add a column that will act as a dummy subset variable
      # and change subset_vars to be the name of this dummy subset variable.
      # This will allow us to use one loop for either case.
      pData(target_data_object)[["DummySubsetVar"]] <- "DummyLevel"
      pData(target_data_object)[["DummySubsetVar"]] <- as.factor(pData(target_data_object)[["DummySubsetVar"]])
      subset_var <- "DummySubsetVar"
    }
    
    # Check that the current subset_var is not the same as either test_var or random_intercept_var.
    if(test_var==subset_var || random_intercept_var==subset_var) {
      print(paste0("Your current subset variable, ", subset_var, " is the same as either your test (contrast) variable or random intercept variable. Please change the subset variable."))
      next
    }
      
    # Get the levels of the current subset_var.
    subset_var_levels <- pData(target_data_object)[[subset_var]] %>% levels
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      # Get all the samples belonging to the current subset_var_level.
      ind <- pData(target_data_object)[[subset_var]] == subset_var_level
      
      # Set the model formula based on whether we have a random slope or not. 
      model_formula <- ifelse(random_slope_i %in% c("no", "FALSE"), "~ TestVar + (1 | RandomInterceptVar)", "~ TestVar + (1 + TestVar | RandomInterceptVar)") %>% as.formula
        
      # Generate the model.
      mixedOutmc <- mixedModelDE(target_data_object[, ind],
                                 elt = "log_norm",
                                 modelFormula = model_formula,
                                 groupVar = "TestVar",
                                 nCores = parallel::detectCores(),
                                 multiCore = FALSE)
      
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
      r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
      r_test <- r_test[, c("Gene", "Subset variable", "Subset level", "Contrast", "Estimate", 
                           "Pr(>|t|)", "FDR")]
      r_test$`Contrast variable` <- test_var
      r_test$`Model number` <- i
      r_test <- r_test %>% dplyr::relocate(`Contrast variable`, .before = Contrast)
      results2 <- rbind(results2, r_test)
      
    }
    
  }
  
  # Reset the variable names. 
  # First check if there's an additional column added if there are no subset variables (see control switch 1a.)
  # If there is, remove it. 
  if("DummySubsetVar" %in% colnames(pData(target_data_object))) pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-DummySubsetVar)
  # Then reset the variable names. 
  colnames(pData(target_data_object)) <- orig_var_names 
  
  # If we changed subset_vars to "DummySubsetVar", change it back to NA.
  if(sum(subset_vars=="DummySubsetVar") == 1) subset_vars <- NA
  
} # End loop level 1: model.


## ---------------------------
# pick top genes for either side of volcano to label
# order genes for convenience:
results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
top_g_list <- list()

model_numbers <- results2$`Model number` %>% unique
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it subset_vars because we already have a variable by that name. ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model_number in model_numbers) {
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model number`==model_number) %>% .$Contrast %>% unique
      for(contrast in contrasts) {
        top_g <- c()
        ind <- results2$`Subset variable`==subset_var & results2$`Subset level`==subset_var_level & results2$`Model number`==model_number & results2$Contrast==contrast
        top_g <- c(top_g,
                   results2[ind, 'Gene'][
                     order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
                   results2[ind, 'Gene'][
                     order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]]) %>%
          unique
        top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]] <- top_g
      }
    }
  }
}

results2 <- results2[, -1*ncol(results2)] # remove invert_P from matrix

## ---------------------------
# Graph

# Categorize results2 based on P-value & FDR for plotting
results2$Color <- "NS or FC < 0.5"
results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
results2$Color <- factor(results2$Color,
                         levels = c("NS or FC < 0.5", "P < 0.05",
                                    "FDR < 0.05", "FDR < 0.001"))
# Set the significance colors.
signif_cols <- c("dodgerblue",
                  "lightblue",
                  "orange2",
                  "grey")
names(signif_cols) <- c("FDR < 0.001", 
                          "FDR < 0.05",  
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
      for(contrast in contrasts) {
        # Get the parameters for this experiment.
        test_var <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Test variable`) %>% unlist %>% .[1]
        random_slope_i <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Random slope`) %>% unlist %>% .[1]
        random_slope_status <- ifelse(random_slope_i %in% c("no", "FALSE"), " | No random slope", " | With random slope")
        test_var_lv_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1] #pData(target_data_object)[[test_var]] %>% levels %>% .[1]
        test_var_lv_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2] #pData(target_data_object)[[test_var]] %>% levels %>% .[2]
        if(subset_var=="NA" || is.na(subset_var) || subset_var=="DummySubsetVar") {
          subset_by <- ""
        } else {
          subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
        }
        
        # Get the top DE genes for this experiment.
        top_genes <- c()
        for(subset_var_level in subset_var_levels) {
          top_genes <- c(top_genes, top_g_list[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]])
        }
        top_genes <- unique(top_genes)
        
        # Plot.
        results2_sub <- results2 %>% dplyr::filter(
          `Subset variable`==subset_var &
            `Subset level`==subset_var_level &
            `Model number`==model_number &
            Contrast==contrast
        )
        plot <- ggplot(results2_sub,
                       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                           color = Color, label = Gene)) +
          geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
          geom_hline(yintercept = -log10(0.05), lty = "dashed") +
          geom_point() +
          scale_color_manual(values = signif_cols) +
          labs(x = paste0("Enriched in ", test_var_lv_1, " <- log2(FC) -> Enriched in ", test_var_lv_2),
               y = "", # Significance, -log10(P)
               # title = paste0("DE genes", subset_by, " \nTest variable: ", test_var, random_slope_status),
               color = "Significance") +
          scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
          geom_text_repel(data = subset(results2_sub, Gene %in% top_genes & FDR < 0.001), # The way we have the graphing for the DEGs set up, it will label all genes for a given subset variable, across all values of that variable. This is because we use the values of the subset variable to facet, and AFAIK, ggplot2 doesn't have a way to exclude labels by the variable that's being faceted by. This may change.
                          size = 4, point.padding = 0.15, color = "black",
                          min.segment.length = .1, box.padding = .2, lwd = 2,
                          max.overlaps = 50) +
          theme_bw(base_size = 16) +
          theme(legend.position = "bottom",
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) # "bottom"
        
        # Add to the list.
        plot_list_diff_exprs[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]] <- plot
      }
    }
  }
}

# Arrange plots into grid.
plot_list_diff_exprs_grid <- list()
for(sv in names(plot_list_diff_exprs)) {
  for(svl in names(plot_list_diff_exprs[[sv]])) {
    for(model_num in names(plot_list_diff_exprs[[sv]][[svl]])) {
      # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
      # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
      # https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0
      
      p_list <- plot_list_diff_exprs[[sv]][[svl]][[model_num]]
      
      n <- length(p_list)
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
      if(subset_var=="NA" || is.na(subset_var) || subset_var=="DummySubsetVar") {
        subset_by <- ""
      } else {
        subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
      }
      
      # Arrange plots in p_list onto a grid.
      plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
      plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Significance, -log10(P)", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = scaling_factor)),
                                                         bottom = grid::textGrob("", gp = grid::gpar(cex = scaling_factor)),
                                                         top = grid::textGrob(paste0("DE genes ", subset_by, " \nTest (contrast) variable: ", test_var, random_slope_status), gp = grid::gpar(cex = scaling_factor)))
      
      # Add back in the legend we extracted earlier. 
      plot_grid2 <- grid.arrange(plot_grid, legend, ncol = 1, heights=c(10, 1))
      
      # Save to list.
      plot_list_diff_exprs_grid[[sv]][[svl]][[model_num]][["plot"]] <- plot_grid2
      plot_list_diff_exprs_grid[[sv]][[svl]][[model_num]][["nCol"]] <- nCol
    }
  }
}

# Add to the PowerPoint.
# Graphing parameters.
plot_width <- 12
plot_height <- 12 
units <- "in"
res <- 300
# Add the graphs.
for(sv in names(plot_list_diff_exprs_grid)) {
  for(svl in names(plot_list_diff_exprs_grid[[sv]])) {
    for(item in names(plot_list_diff_exprs_grid[[sv]][[svl]])) {
      plot <- plot_list_diff_exprs_grid[[sv]][[svl]][[item]][["plot"]] %>% ggpubr::as_ggplot()
      nCol <- plot_list_diff_exprs_grid[[sv]][[svl]][[item]][["nCol"]]
      
      # Set the scaling factors.
      # It's the same for both width and height since the grids are squares.
      scaling_factor <- ifelse(nCol > 1, (nCol / 2)^2, 1)  # Number of rows in current grid / 2 (base number)
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, "LMM-differential-expression_graph-", sv, "-", svl, "-", item, ".eps")
      png_path <- paste0(output_dir_pubs, "LMM-differential-expression_graph-", sv, "-", svl, "-", item, ".png")
      saveEPS(plot, eps_path, width = (plot_width * scaling_factor), height = (plot_height * scaling_factor))
      savePNG(plot, png_path, width = (plot_width * scaling_factor), height = (plot_height * scaling_factor), units = units, res = (res * scaling_factor))

      # Add to the PowerPoint.
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Differential expression"),
                         location = ph_location_label(ph_label = "Title 1")) %>%
        officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                         location = ph_location_label(ph_label = "Content Placeholder 2"),
                         use_loc_size = FALSE)
    }
  }
}

## ---------------------------
# Export to disk.

# Export NanoStringGeoMxSet.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_differential-expression.rds"))
# Export table of DE genes to CSV.
results2 %>% write.csv(paste0(output_dir_tabular, "LMM-differential-expression_results.csv")) 
# Output everything to the PowerPoint. 
print(pptx, cl_args[5])