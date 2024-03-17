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

# Convert test variables to factors.
for(test_var in test_vars) {
  pData(target_data_object)[[test_var]] <- pData(target_data_object)[[test_var]] %>% as.factor
}

# Create the table of LMM parameter combinations.
param_combos <- cbind(random_slope, test_vars, random_intercept_vars, random_slope_vars) %>% as.data.frame %>% dplyr::mutate(`Model number` = 1:nrow(.))
colnames(param_combos)[1:4] <- c("Random slope", "Test variable", "Random intercept variable", "Random slope variable")
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
# Now we can loop over all the param combos and run LMM on each one.
for(i in 1:nrow(param_combos)) {
  # Get the params for this experiment. 
  random_slope_i <- param_combos[i,1]
  test_var <- param_combos[i,2]
  random_intercept_var <- param_combos[i,3]
  random_slope_vars <- param_combos[i,4]

  # Reset the variable names.
  colnames(pData(target_data_object)) <- orig_var_names
  
  # Rename the test variable (test_var) and the random intercept variable (random_intercept_var).
  pData(target_data_object) <- pData(target_data_object) %>% 
    dplyr::rename(TestVar = !!as.name(test_var), RandomInterceptVar = !!as.name(random_intercept_var))
  
  if(random_slope_i %in% c("no", "FALSE")) {
    mixedOutmc <- mixedModelDE(target_data_object,
                                      elt = "log_norm",
                                      modelFormula = ~ TestVar + (1 | RandomInterceptVar),
                                      groupVar = "TestVar",
                                      nCores = parallel::detectCores(),
                                      multiCore = FALSE)
                         
    
  } else {
    mixedOutmc <- mixedModelDE(target_data_object,
                 elt = "log_norms",
                 modelFormula = ~ TestVar + (1 + TestVar | RandomInterceptVar),
                 groupVar = "TestVar",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
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
  # r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")] # r_test[, c("Gene", "Subset", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
  r_test$`Contrast variable` <- test_var
  r_test$`Model number` <- i
  r_test <- r_test %>% dplyr::relocate(`Contrast variable`, .before = Contrast)
  results2 <- rbind(results2, r_test)
  
  # Reset the variable names.
  colnames(pData(target_data_object)) <- orig_var_names 

}

## ---------------------------
# Graph

# Initialize the list to hold the plots.
plot_list_diff_exprs <- list()

# Categorize results2 based on P-value & FDR for plotting
results2$Color <- "NS or FC < 0.5"
results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
results2$Color <- factor(results2$Color,
                         levels = c("NS or FC < 0.5", "P < 0.05",
                                    "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
top_g_list <- list()
# for(cond in unique(pData(target_data_object)[[grouping_var]])) { # grouping_var not yet defined 2024/03/11.
#     ind <- results2$Subset == cond
#     top_g <- c(top_g,
#                results2[ind, 'Gene'][
#                    order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
#                results2[ind, 'Gene'][
#                    order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]])
# }
for(model_number in unique(results2$`Model number`)) {
  top_g <- c()
  ind <- results2$`Model number` == model_number
  top_g <- c(top_g,
             results2[ind, 'Gene'][
               order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             results2[ind, 'Gene'][
               order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]]) %>%
    unique
  top_g_list[[test_var]] <- top_g
}
results2 <- results2[, -1*ncol(results2)] # remove invert_P from matrix

# Graph results2
for(model_number in unique(results2$`Model number`)) {
  test_var <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Test variable`) %>% unlist %>% .[1]
  random_slope_i <- param_combos %>% dplyr::filter(`Model number`==model_number) %>% dplyr::select(`Random slope`) %>% unlist %>% .[1]
  random_slope_status <- ifelse(random_slope_i %in% c("no", "FALSE"), " - no random slope", " - random slope")

  test_var_lv_1 <- pData(target_data_object)[[test_var]] %>% levels %>% .[1]
  test_var_lv_2 <- pData(target_data_object)[[test_var]] %>% levels %>% .[2]
  plot_list_diff_exprs[[paste0("model_", model_number)]] <- ggplot(results2,
                                             aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                                 color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = paste0("Enriched in ", test_var_lv_1, " <- log2(FC) -> Enriched in ", test_var_lv_2),
         y = "Significance, -log10(P)",
         title = paste0("DE genes - ", test_var, random_slope_status),
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results2, Gene %in% top_g & FDR < 0.001),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") #+
  # facet_wrap(~Subset, scales = "free_y")
}

# Export NanoStringGeoMxSet.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_differential-expression.rds"))
# Export table of DE genes to CSV.
results2 %>% write.csv(paste0(output_dir_tabular, "LMM-differential-expression_results.csv")) 

# Add to the PowerPoint.
# Add the graphs.
for(item in names(plot_list_diff_exprs)) {
  # Add to the PowerPoint. 
  pptx <- pptx %>%
    officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
    officer::ph_with(value = paste0("Differential expression"),
                      location = ph_location_label(ph_label = "Title 1")) %>% 
    officer::ph_with(value = plot_list_diff_exprs[[item]],
                      location = ph_location_label(ph_label = "Content Placeholder 2"))
}
print(pptx, cl_args[5])