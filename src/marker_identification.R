## Source the setup.R file.
source("src/setup.R")

###################################################################
##
## Setup 
##
###################################################################

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

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

# Add to the PowerPoint. 
# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Marker identification (linear mixed models)"), 
                   location = ph_location_label(ph_label = "Title 1"))

###################################################################
##
## Heatmaps of differentially expressed genes
##
###################################################################

## ----------------------------------------------------------------
##
## Run LMM
##
## ----------------------------------------------------------------
# Formula follows conventions defined by the lme4 package.
# When running LMM, mixedModelDE seems to have issues with variable names with spaces etc., even if enclosed in backticks (``). 
# So we're going to rename the variables of interest to something without spaces and other special characters.
# But first, we'll have to save the original variable names so every time we rename the variables in the for() loop,
# we can reset them afterwards.
orig_var_names <- colnames(pData(target_data_object))
results3 <- c()
markers <- list()

# Loop levels (outermost to innermost):
# 1) Model (param_combos)
# 2) Subset variable (subset_vars)

# Now we can loop over all the param combos and run LMM on each one.
# Loop level 1: model.
for(i in 1:nrow(param_combos)) {
  markers[[paste0("model_", i)]] <- list()
  
  # Make a copy of the original target_data_object.
  target_data_object_2 <- target_data_object
  
  print(paste0("Working on model #", i, "."))
  
  # Get the params for this experiment. 
  random_slope_i <- param_combos[i,1]
  test_var <- param_combos[i,2]
  random_intercept_var <- param_combos[i,3]
  # random_slope_vars <- param_combos[i,4]
  
  # Reset the variable names.
  colnames(pData(target_data_object_2)) <- orig_var_names
  
  # Rename the test variable (test_var) and the random intercept variable (random_intercept_var).
  pData(target_data_object_2) <- pData(target_data_object_2) %>%
    dplyr::rename(TestVar = !!as.name(test_var), RandomInterceptVar = !!as.name(random_intercept_var))
  
  # control switch 1 (subset variables exist or not) << loop level 1 (model).
  if(sum(is.na(subset_vars)) == length(subset_vars) || sum(subset_vars=="NA", na.rm = T) == length(subset_vars[!is.na(subset_vars)])) {
    # control switch 1a: no subset variables << loop level 1 (model).
    
    # Since there are no subset variables, we will add a column that will act as a dummy subset variable
    # and change subset_vars to be the name of this dummy subset variable.
    # This will allow us to use one loop for either case (controls switch 1a or 1b).
    pData(target_data_object_2)[["DummySubsetVar"]] <- "DummyLevel"
    pData(target_data_object_2)[["DummySubsetVar"]] <- as.factor(pData(target_data_object_2)[["DummySubsetVar"]])
    subset_vars <- c("DummySubsetVar")
    
  } # End control switch 1a (no subset variables) << loop level 1 (model).
  
  # loop level 2 (subset variable) << loop level 1 (model)
  for(subset_var in subset_vars) {
    print(paste0("Working on subset variable ", subset_var, " for model #", i, "."))
    
    # Check if the current subset_var is NA.
    if(subset_var == "NA" || is.na(subset_var)) {
      # Add a column that will act as a dummy subset variable
      # and change subset_vars to be the name of this dummy subset variable.
      # This will allow us to use one loop for either case.
      pData(target_data_object_2)[["DummySubsetVar"]] <- "DummyLevel"
      pData(target_data_object_2)[["DummySubsetVar"]] <- as.factor(pData(target_data_object_2)[["DummySubsetVar"]])
      subset_var <- "DummySubsetVar"
      markers[[paste0("model_", i)]][[subset_var]] <- list()
    }
    
    # Check that the current subset_var is not the same as either test_var or random_intercept_var.
    if(test_var==subset_var || random_intercept_var==subset_var) {
      print(paste0("The current subset variable, ", subset_var, " is the same as either your test (contrast) variable or random intercept variable. Skipping this subset variable - contrast/random intercept variable combination."))
      next
    }
    
    markers[[paste0("model_", i)]][[subset_var]] <- list() # Putting it here otherwise "NA" will be one of the names of the list items. 
    # Get the levels of the current subset_var.
    subset_var_levels <- pData(target_data_object_2)[[subset_var]] %>% as.factor %>% levels # as.factor needed because it might be a character vector.
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      markers[[paste0("model_", i)]][[subset_var]][[subset_var_level]] <- list()
      print(paste0("Working on level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
      
      # Get all the samples belonging to the current subset_var_level.
      ind <- pData(target_data_object_2)[[subset_var]] == subset_var_level
      
      # Set the model formula based on whether we have a random slope or not. 
      model_formula <- ifelse(random_slope_i %in% c("no", "FALSE"), "~ TestVar + (1 | RandomInterceptVar)", "~ TestVar + (1 + TestVar | RandomInterceptVar)") %>% as.formula
      
      # loop level 4 (current normalization method) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
      for(norm_method in normalization_methods) {
        markers[[paste0("model_", i)]][[subset_var]][[subset_var_level]][[norm_method]] <- list()
        print(paste0("Working on normalization method ", norm_method, " for level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
        
        # Create a log2 transform of the data for analysis.
        assayDataElement(object = target_data_object_2, elt = "log_norm", validate = FALSE) <- # Have to set validate to FALSE; otherwise it thinks the dimensions aren't the same. ... 
          assayDataApply(target_data_object_2, 2, FUN = function(x) log2(x+1), elt = norm_method)
        
        # loop level 5 (level of current test variable) << loop level 4 (current normalization method) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
        test_var_levels <- pData(target_data_object_2) %>% .$TestVar %>% levels
        for(test_var_level in test_var_levels) {
          print(paste0("Working on level ", test_var_level, " for the test variable; normalization method ", norm_method, " for level ", subset_var_level, " for subset variable ", subset_var, " for model #", i, "."))
          
          # Set the random seed.
          set.seed(random_seed)
          
          # Create a temporary target data object in which
          # the levels of the test variable that are NOT the current test_var_level
          # will be set to "Other"; the comparison will then be
          # the current test_var_level vs. "Other".
          tdo_tmp <- target_data_object_2 %>% subset(!(TargetName %in% neg_probes), ind) # Exclude negative probes (b/c variance might be 0, and this will screw up LLM); and include only the samples in ind (samples belonging to current subset_var_level, loop level 3.)
          pData(tdo_tmp)$TestVar <- as.character(pData(tdo_tmp)$TestVar)
          pData(tdo_tmp)$TestVar[pData(tdo_tmp)$TestVar != test_var_level] <- "ZZZ" # Make sure the other levels aren't alphabetically after this! 
          pData(tdo_tmp)$TestVar <- as.factor(pData(tdo_tmp)$TestVar)
          
          # Generate the model.
          mixedOutmc <- mixedModelDE(tdo_tmp,
                                     elt = "log_norm",
                                     modelFormula = model_formula,
                                     groupVar = "TestVar",
                                     nCores = parallel::detectCores(),
                                     multiCore = TRUE)
          rm(tdo_tmp)
          
          # Format results as data.frame.
          r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
          tests <- rownames(r_test)
          r_test <- as.data.frame(r_test)
          r_test$Contrast <- tests
          r_test$`Contrast level` <- test_var_level
          
          # Use lapply in case you have multiple levels of your test factor to
          # correctly associate gene name with its row in the results table.
          r_test$Gene <- 
            unlist(lapply(colnames(mixedOutmc),
                          rep, nrow(mixedOutmc["lsmeans", ][[1]])))
          r_test$`Subset variable` <- subset_var
          r_test$`Subset level` <- subset_var_level
          r_test$`Normalization method` <- normalization_names[names(normalization_names)==norm_method]
          r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
          r_test <- r_test[, c("Gene", "Subset variable", "Subset level", "Normalization method", "Contrast", "Contrast level", "Estimate", 
                               "Pr(>|t|)", "FDR")]
          r_test$`Contrast variable` <- test_var
          r_test$`Model number` <- i
          r_test <- r_test %>% dplyr::relocate(`Contrast variable`, .before = Contrast)
          results3 <- rbind(results3, r_test)
          
          # Add the markers to the vector. 
          markers[[paste0("model_", i)]][[subset_var]][[subset_var_level]][[norm_method]][[test_var_level]] <- 
            r_test %>% dplyr::filter(FDR < de_genes_cutoffs[1]) %>% .$Gene
          
        } # End loop level 5: test variable levels.
        
      } # End loop level 4: normalization method. 
      
    }
    
  }
  
  # Reset the variable names. 
  # First check if there's an additional column added if there are no subset variables (see control switch 1a.)
  # If there is, remove it. 
  if("DummySubsetVar" %in% colnames(pData(target_data_object))) pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-DummySubsetVar)
  # Then reset the variable names. 
  colnames(pData(target_data_object)) <- orig_var_names 
  
  # If we changed subset_vars to "DummySubsetVar", change it back to NA.
  if(sum(subset_vars=="DummySubsetVar", na.rm = T) == 1) subset_vars <- NA
  
  # Remove target data object copy.
  rm(target_data_object_2)
  gc()
  
} # End loop level 1: model.


## ----------------------------------------------------------------
##
## Create heatmaps
##
## ----------------------------------------------------------------
de_heatmaps <- list()
for(model_num in names(markers)) {
  de_heatmaps[[model_num]] <- list()
  i <- model_num %>% regexPipes::gsub("model_", "") %>% as.integer
  test_var <- param_combos[i,2]
  
  for(subset_var in names(markers[[model_num]])) {
    if(length(markers[[model_num]][[subset_var]]) < 1) next
    if(subset_var %in% c("NA", "DummySubsetVar")) next
    
    de_heatmaps[[model_num]][[subset_var]] <- list()
    
    for(subset_var_level in names(markers[[model_num]][[subset_var]])) {
      de_heatmaps[[model_num]][[subset_var]][[subset_var_level]] <- list()
      # Get all the samples belonging to the current subset_var_level.
      if(subset_var == "DummySubsetVar") {
        ind <- 1:nrow(pData(target_data_object))
      } else {
        ind <- pData(target_data_object)[[subset_var]] == subset_var_level
      }
      
      
      for(norm_method in names(markers[[model_num]][[subset_var]][[subset_var_level]])) {
        print(paste0("Creating heatmap for model number ", i, ", subset variable ", subset_var, ", subset variable level ", subset_var_level, ", normalization method ", norm_method))
        
        # Make a copy of the original target_data_object.
        target_data_object_2 <- target_data_object
        
        # Create a log2 transform of the data for analysis.
        assayDataElement(object = target_data_object_2, elt = "log_norm", validate = FALSE) <- # Have to set validate to FALSE; otherwise it thinks the dimensions aren't the same. ... 
          assayDataApply(target_data_object_2, 2, FUN = function(x) log2(x+1), elt = norm_method)
        
        GOI <- markers[[model_num]][[subset_var]][[subset_var_level]][[norm_method]] %>% 
          unlist() %>% 
          unique() %>% 
          .[!(. %in% neg_probes)]
        exprs_mat <- assayDataElement(target_data_object_2[GOI, ind], elt = "log_norm")
        annot <- data.frame(pData(target_data_object_2)[ind, heatmap_ann_vars])
        rownames(annot) <- colnames(exprs_mat)
        de_heatmaps[[model_num]][[subset_var]][[subset_var_level]][[norm_method]][["nMarkers"]] <- length(GOI)
        de_heatmaps[[model_num]][[subset_var]][[subset_var_level]][[norm_method]][["heatmap"]] <- pheatmap(exprs_mat,
                                                                                                           scale = "row",
                                                                                                           show_rownames = TRUE, show_colnames = FALSE,
                                                                                                           border_color = NA,
                                                                                                           clustering_method = "average",
                                                                                                           clustering_distance_rows = "correlation",
                                                                                                           clustering_distance_cols = "correlation",
                                                                                                           breaks = seq(-3, 3, 0.05),
                                                                                                           color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
                                                                                                           annotation_col = annot, # https://www.researchgate.net/post/R-error-gpar-element-fill-must-not-be-length-0
                                                                                                           # Apparently pData loses the rownames, hence using the colnames of exprs_mat as the rownames of the annotation data frame.
                                                                                                           main = paste0("Test variable: ", test_var, 
                                                                                                                         "\nSubset variable: ", subset_var, " | Subset level: ", subset_var_level,
                                                                                                                         "\nNormalization method: ", norm_method)
        )
        
        rm(target_data_object_2)
        gc()
      }
    }
  }
}

## ----------------------------------------------------------------
##
## Arrange heatmaps into grid
##
## ----------------------------------------------------------------
# plot_list_de_heatmap_grids <- list()
# for(model_num in names(de_heatmaps)) {
# 
#   for(subset_var in names(de_heatmaps[[model_num]])) {
#     
#     for(subset_var_level in names(de_heatmaps[[model_num]][[subset_var]])) {
#       
#       p_list <- de_heatmaps[[model_num]][[subset_var]][[subset_var_level]]
#       
#       n <- length(p_list)
#       nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
#       
#       # Set the scaling factors for label and legend size.
#       sqrt_n_col <- sqrt(nCol)
#       scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 2), 1) # Number of rows in current grid / 2 (base number)
#       
#       # Get the model information (test/contrast variable, random slope status, etc.)
#       i <- model_num %>% regexPipes::gsub("model_", "")
#       random_slope_i <- param_combos[i,1]
#       test_var <- param_combos[i,2]
#       random_intercept_var <- param_combos[i,3]
#       random_slope_status <- ifelse(random_slope_i %in% c("no", "FALSE"), " | No random slope", " | With random slope")
#       if(subset_var=="NA" || is.na(subset_var) || subset_var=="DummySubsetVar") {
#         subset_by <- ""
#       } else {
#         subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
#       }
#       
#       # Arrange plots in p_list onto a grid.
#       plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
#       plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Significance, -log10(P)", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = scaling_factor)),
#                                                          bottom = grid::textGrob("", gp = grid::gpar(cex = scaling_factor)),
#                                                          top = grid::textGrob(paste0("DE genes ", subset_by, 
#                                                                                      " \nTest (contrast) variable: ", test_var, random_slope_status), 
#                                                                               gp = grid::gpar(cex = scaling_factor)))
#       
#       # Add back in the legend we extracted earlier. 
#       plot_grid2 <- grid.arrange(plot_grid, legend, ncol = 1, heights=c(10, 1))
#       
#       # Save to list.
#       plot_list_diff_exprs_grid[[subset_var]][[subset_var_level]][[model_num]][[normalization_method]][["plot"]] <- plot_grid2
#       plot_list_diff_exprs_grid[[subset_var]][[subset_var_level]][[model_num]][[normalization_method]][["nCol"]] <- nCol
#     
#     }
#   }
# }

## ----------------------------------------------------------------
##
## Save to EPS/PNG and then write to PowerPoint
##
## ----------------------------------------------------------------
# Graphing parameters.
plot_width <- 12
plot_height <- 12 
units <- "in"
res <- 300
for(model_num in names(de_heatmaps)) {
  i <- model_num %>% regexPipes::gsub("model_", "") %>% as.integer
  test_var <- param_combos[i,2]
  
  for(subset_var in names(de_heatmaps[[model_num]])) {
    
    for(subset_var_level in names(de_heatmaps[[model_num]][[subset_var]])) {
      
      for(norm_method in names(de_heatmaps[[model_num]][[subset_var]][[subset_var_level]])) {
        plot <- de_heatmaps[[model_num]][[subset_var]][[subset_var_level]][[norm_method]][["heatmap"]]
        nMarkers <- de_heatmaps[[model_num]][[subset_var]][[subset_var_level]][[norm_method]][["nMarkers"]]
        
        # Set the scaling factor--it's used only for height.
        scaling_factor_y <- (nMarkers / 50)  # Number of rows in current grid / 2 (base number)
        addition_factor_y <- ifelse(nMarkers < 25, 2, 0) # Add 2 inches to cover the non-heatmap parts so scaling_factor_y doesn't make the heatmap too small.
        # scaling_factor_x <- ifelse(nMarkers <= 25, nMarkers / 50, 1)
        # scaling_factor_res <- mean(c(scaling_factor_y, scaling_factor_x))
        
        # Save to EPS and PNG and then ...
        eps_path <- paste0(output_dir_pubs, "LMM-marker-genes_heatmap-", test_var, "-", subset_var, "-", subset_var_level, "-", norm_method, ".eps")
        png_path <- paste0(output_dir_pubs, "LMM-marker-genes_heatmap-", test_var, "-", subset_var, "-", subset_var_level, "-", norm_method, ".png")
        saveEPS(plot, eps_path, width = plot_width, height = ((plot_height * scaling_factor_y) + addition_factor_y))
        savePNG(plot, png_path, width = plot_width, height = ((plot_height * scaling_factor_y) + addition_factor_y), units = units, res = (res))
        
        # Add to the PowerPoint.
        pptx <- pptx %>%
          officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
          officer::ph_with(value = paste0("Marker genes"),
                           location = ph_location_label(ph_label = "Title 1")) %>%
          officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                           location = ph_location_label(ph_label = "Content Placeholder 2"),
                           use_loc_size = FALSE)
        
      }
    }
  }
}

###################################################################
##
## Write everything to disk.
##
###################################################################

# Export NanoStringGeoMxSet.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_marker-identification.rds"))
# Export graphs.
saveRDS(de_heatmaps, paste0(output_dir_rdata, "LMM-marker_heatmaps.rds"))
# Export tables of DE genes to CSV.
results3 %>% write.csv(paste0(output_dir_tabular, "LMM-marker_results.csv"))
# Output everything to the PowerPoint. 
print(pptx, cl_args[5])