## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

# Set the normalization method
normalization_method <- normalization_names[names(normalization_names)==normalization_methods[1]]

# Cell profile matrices for SpatialDecon
data("safeTME")
data("safeTME.matches")

# # Set path to CIBERSORT required files
# set_cibersort_binary(path_to_cibersort)
# set_cibersort_mat(path_to_lm22)
# 
# # Calculate TPM - this is necessary for CIBERSORT among others, not so much for xCell or MCP-counter.
# # https://bioinformatics.stackexchange.com/questions/2567/how-can-i-calculate-gene-length-for-rpkm-calculation-from-counts-data
# # But can we even calculate TPM for GeoMx data? Bc it's probe-based, so it wouldn't have the same assumptions that RNA-seq does ...

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Immune deconvolution ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## ................................................
##
### Calculation ----
##
## ................................................
# Subset the target data object to include only the transcriptome probes.
target_data_object_exprs <- subset(target_data_object, Module %in% main_module)
# Extract the expression matrix.
exprs_mat <- target_data_object_exprs@assayData[[normalization_methods[1]]]
# If the species is mouse (Mus musculus), we can still run human methods on mouse data, but we just need to convert to their human orthologues.
# Also, download the mouse profile matrix using SpatialDecon::download_profile_matrix()
if(species == "Mus musculus") {
  exprs_mat_ortho <- convert_human_mouse_genes(exprs_mat, convert_to = 'human') # convert_human_mouse_genes() in src/helper_functions.R
  system.time(download_profile_matrix(species = "Mouse", age_group = "Adult", matrixname = "ImmuneAtlas_ImmGen")) # https://github.com/Nanostring-Biostats/CellProfileLibrary/blob/master/Mouse/Mouse_datasets_metadata.csv
  # user  system elapsed 
  # 0.138   0.059   1.858 
}

# Loop through imm_decon_methods.
imm_decon_res_list <- list()
for(method in imm_decon_methods) {
  if(!(method %in% c("quantiseq", "mcp_counter", "xcell", "epic", "abis", "estimate", "spatialdecon", "mmcp_counter"))) {
    warning(glue::glue("{method} is currently not supported by this pipeline. Please note that TIMER and ConsensusTME, while included in the immunedeconv package, are currently not available in this pipeline due to extra arguments that must be passed to the function; and CIBERSORT will not be available until we figure out how to make the source code play nicely with the immunedeconv package."))
    next
  } else {
    # Error catching: https://stackoverflow.com/a/55937737/23532435
    skip_to_next <- FALSE
    
    # If the method is mmcp_counter, check if the species is mouse (Mus musculus).
    if(method == "mmcp_counter") {
      if(species != "Mus musculus") {
        warning("mMCP-counter can only be used with mouse data. Skipping to the next one.")
        next
      } else {
        # Using a mouse method on mouse data. Do not convert to orthologues.
        exprs_mat_effective <- exprs_mat
      }
    } else {
      # If the current method is not mMCP-counter, check the species.
      if(species == "Mus musculus") {
        # Using a human method on mouse data. Convert to orthologues.
        exprs_mat_effective <- exprs_mat_ortho
      } else {
        # Using a human method on human data. Do not convert to orthologues.
        exprs_mat_effective <- exprs_mat
      }
    }
    
    # Special steps needed for SpatialDecon
    if(method == "spatialdecon") {
      # The spatialdecon function takes 3 arguments of expression data:
      #   
      # 1. The normalized data.
      # 2. A matrix of expected background for all data points in the normalized data matrix.
      # 3. Optionally, either a matrix of per-data-point weights, or the raw data, which is used to derive weights (low counts are less statistically stable, and this allows spatialdecon to down-weight them.)
      # https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html
      
      # If the user provided a custom profile matrix, load it
      load_safeTME <- TRUE
      if(!flagVariable(spatial_decon_profile_matrix)) {
        if(file.exists(spatial_decon_profile_matrix)) {
          load_spatial_decon_profile_matrix <- tryCatch(
                                                        {load(spatial_decon_profile_matrix)}, 
                                                        error = function(msg) {return(NA)}
                                                        )
          if(exists("profile_matrix")) load_safeTME <- FALSE
        }
      }
      if(load_safeTME) { 
        message(glue::glue("The default reference matrix, safeTME, will be used with the SpatialDecon method. If this was not your intended reference matrix, you may have either 1) left the `spatial_decon_profile_matrix` field in the configuration YAML file blank or 2) did not fill in a full path to an .RData file containing an object named `profile_matrix`: a matrix with genes as rows and cell types as columns. If you did intend to use safeTME, you can safely ignore this message"))
        profile_matrix <- safeTME 
      }
      
      # Estimate each data point's expected BG from the negative control probes from its corresponding observation.
      # Also, the SpatialDecon algorithm will do its own background subtraction
      # so we have to give it data that has been normalized using a non BG-sub method (otherwise we'll get a LinAlg error)
      # we'll go with quant
      negnames <- fData(target_data_object_exprs) %>% 
        dplyr::filter(Negative == TRUE & Module %in% modules_exprs) %>% 
        .$TargetName
      bg <- derive_GeoMx_background(norm = target_data_object_exprs@assayData[["quant"]],
                                    probepool = fData(target_data_object_exprs)$Module,
                                    negnames = negnames)
      
      signif(profile_matrix[seq_len(3), seq_len(3)], 2)
      # heatmap(sweep(safeTME, 1, apply(safeTME, 1, max), "/"),
      #         labRow = NA, margins = c(10, 5))
      
      # Set the cell profile matrix based on the species.
      if(species == "Homo sapiens" & !exists("profile_matrix")) {
        cpm <- safeTME
        message("Using safeTME as reference (profile) matrix for SpatialDecon")
      } else {
        cpm <- profile_matrix %>% as.matrix
        message("Using custom reference (profile) matrix for SpatialDecon")
      }
      
      # Run spatial deconvolution.
      system.time({
        res <- tryCatch(runspatialdecon(object = target_data_object_exprs,
                                        norm_elt = "quant", # "neg_norm"   "log_norm"   "bg_sub"     "exprs"      "bg_sub_neg" "quant"      "bg_sub_q3"  "q3_norm"   
                                        raw_elt = "exprs",
                                        X = cpm,
                                        align_genes = TRUE#,
                                        #cell_counts = pData(target_data_object_exprs)
                                        ),
                        error = function(e) {skip_to_next <<- TRUE})
      })
      if(class(res) != "logical") imm_decon_res <- res$beta %>% t %>% as.data.frame %>% rownames_to_column("cell_type")
    } else {
      imm_decon_res <- tryCatch(immunedeconv::deconvolute(exprs_mat_effective, method),
                                error = function(e) {skip_to_next <<- TRUE})
    }
    
    if(skip_to_next) {
      warning(glue::glue("An error occurred when trying to run immune deconvolution method {method}. Skipping to the next method"))
      next
    }
    
    imm_decon_res_list[[method]] <- imm_decon_res
  }
}

## ................................................
##
### Differential abundance ----
##
## ................................................
if(!flagVariable(lmm_formulae_immune)) {
  # By "differential abundance," we mean _between-sample_ comparisons of immune-cell populations,
  # NOT between-cell-type comparisons.
  # Per https://omnideconv.org/immunedeconv/articles/immunedeconv.html,
  # the following methods allow between-sample comparisons:
  # MCP-counter, xCell, TIMER, ConsensusTME, ESTIMATE, ABIS, mMCP-counter, BASE, EPIC, quanTIseq, CIBERSORT abs., seqImmuCC
  between_sample_methods <- c("mcp_counter", "xcell", "estimate", "abis", "mmcp_counter", "epic", "quantiseq", "spatialdecon")
  model_number <- 1
  da_res_list <- list()
  if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete dataset"
  
  for(formula in lmm_formulae_immune) {
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
    df_cols <- sapply(pData(target_data_object), is.data.frame)
    
    # Convert only non-data-frame columns (except Sample) to factors
    pData_tmp <- pData(target_data_object) %>%
      tibble::rownames_to_column(var = "Sample") %>% 
      mutate(across(!is.data.frame, ~ if (is.numeric(.)) . else as.factor(.))) 
    # %>% select(all_of(formula_vars))
    # Add "Complete dataset" as variable if it doesn't exist already
    if(subset_vars_imm_decon=="Complete dataset" & !("Complete dataset" %in% colnames(pData_tmp))) pData_tmp$`Complete dataset` <- as.factor("Complete dataset")
    
    # Loop level 2 (deconvolution method)
    for(method in names(imm_decon_res_list)) {
      da_res_list[[paste0("model_", model_number)]][[method]] <- list()
      
      # Check if the method can be used in between-sample comparisons
      if(!(method %in% between_sample_methods)) next
      
      # Convert immune_data to long format
      # and merge with `pData_subset`
      immune_long <- imm_decon_res_list[[method]] %>%
        pivot_longer(cols = -cell_type, names_to = "Sample", values_to = "score") %>%
        left_join(pData_tmp, by = "Sample")
      
      # Loop level 3 (subset variable)
      for(subset_var in subset_vars_imm_decon) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete dataset"` above. 
        da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]] <- list()
        
        if(flagVariable(subset_var)) {
          subset_tag <- "All ROIs"
          subset_var <- "Complete dataset"
        } else {
          subset_tag <- subset_var
        }
        
        # Get the levels of the current subset_var.
        subset_var_levels <- pData_tmp[[subset_var]] %>% levels # Previously as.factor %>% needed because it might be a character vector, but now we convert all non-data-frame columns (except Sample) to factor above
        
        # If subset_var_imm_decon_levels_manual is set, filter subset_var_levels to include only those values
        subset_var_imm_decon_levels_manual_i <- subset_var_imm_decon_levels_manual[[subset_var]]
        if(sum(is.na(subset_var_imm_decon_levels_manual_i)) < length(subset_var_imm_decon_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
          if(sum(subset_var_imm_decon_levels_manual_i %in% subset_var_levels) < 1) {
            warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
          } else { # At least one subset_var_level_manual_i is an actual level of the current subset variable
            subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_imm_decon_levels_manual_i]
          }
        }
        
        # loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
        for(subset_var_level in subset_var_levels) {
          da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]][[subset_var_level]] <- list()
          
          # Get all the samples belonging to the current subset_var_level.
          samples <- pData_tmp %>% dplyr::filter(!!as.name(subset_var)==subset_var_level) %>% dplyr::select(Sample) %>% unlist
          # ind <- pData_tmp[[subset_var]] == subset_var_level
          # pData_tmp_sub <- pData_tmp[ind,]
          
          # loop level 5 (cell type) << loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
          for(cell in unique(immune_long$cell_type)) {
            message("")
            message(glue::glue("Calculating differential abundance using formula {formula} | {method} method | subset variable {subset_var} | level {subset_var_level} for cell type {cell} \n"))
            skip_to_next <- FALSE
            
            # Subset data for this cell type, subset variable, and subset variable level
            cell_data <- immune_long %>%
              dplyr::filter((Sample %in% samples) & (cell_type == cell))
            
            # Fit the user-defined linear mixed model
            model <- tryCatch(
              lmerTest::lmer(as.formula(formula_af), data = cell_data),
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
            
          } # End loop level 5 (cell type) << loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
          
        } # End loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
        
      } # End loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
      
      # # Combine results into a data frame
      # results_df <- bind_rows(da_res_list[[paste0("model_", model_number)]][[method]])
      
    } # End loop level 2 (deconvolution methods)
    
    model_number <- model_number + 1
    
  } # End loop level 1 (LMM model)
  
  da_res_df <- bind_rows(rlang::squash(da_res_list)) # `squash` recursively flattens the list.
  
} # End if(!flagVariable(lmm_formulae_immune)) 

## ................................................
##
### Visualization ----
##
## ................................................
# https://omnideconv.org/immunedeconv/articles/detailed_example.html
# Rename the samples to reflect their segment and type.
pData_tmp <- pData(target_data_object) %>% as.data.frame %>% tibble::rownames_to_column(var = "rownames")
observation_identifiers <- intersect(observation_identifiers, colnames(pData_tmp)) # Make sure the observation identifiers are actually in the pData.
pData_tmp <- pData_tmp %>% tidyr::unite("All ROIs", c(observation_identifiers, rownames), remove = FALSE, na.rm = FALSE, sep = " | ") 
pData_tmp$`All ROIs` <- pData_tmp$`All ROIs` %>% regexPipes::gsub("\\.dcc", "")
pData_tmp$`Complete dataset` <- "Complete dataset"
plot_list <- list()
if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete dataset"
for(method in names(imm_decon_res_list)) {
  plot_list[[method]] <- list()
  df <- imm_decon_res_list[[method]]

  # QuanTIseq, CIBERSORT (absolute), Epic - visualize as stacked bar charts.
  # MCP-counter, SpatialDecon - visualize as box/dot/violin plot.
  
  for(subset_var in subset_vars_imm_decon) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete dataset"` above. 
    plot_list[[method]][[subset_var]] <- list()
    
    if(flagVariable(subset_var)) {
      subset_tag <- "All ROIs"
      subset_var <- "Complete dataset"
    } else {
      subset_tag <- subset_var
    }
    
    # Get the levels of the current subset_var.
    subset_var_levels <- pData_tmp[[subset_var]] %>% as.factor %>% levels # as.factor needed because it might be a character vector.
    
    # If subset_var_imm_decon_levels_manual is set, filter subset_var_levels to include only those values
    subset_var_imm_decon_levels_manual_i <- subset_var_imm_decon_levels_manual[[subset_var]]
    if(sum(is.na(subset_var_imm_decon_levels_manual_i)) < length(subset_var_imm_decon_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
      if(sum(subset_var_imm_decon_levels_manual_i %in% subset_var_levels) < 1) {
        warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
      } else { # At least one subset_var_level_manual_i is an actual level of the current subset variable
        subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_imm_decon_levels_manual_i]
      }
    }
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      plot_list[[method]][[subset_var]][[subset_var_level]] <- list()
      
      # Get all the samples belonging to the current subset_var_level.
      ind <- pData_tmp[[subset_var]] == subset_var_level
      pData_tmp_sub <- pData_tmp[ind,]
      
      for(grouping_var in grouping_vars_imm_decon) {
        plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]] <- list()
        df2 <- df
        
        # Gather the dataframe for graphing.
        if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
          df3 <- df2 %>% 
            gather(`All ROIs`, fraction, -cell_type)
        } else {
          df3 <- df2 %>% 
            gather(`All ROIs`, score, -cell_type)
        }
        df3$`All ROIs` <- df3$`All ROIs` %>% regexPipes::gsub("\\.dcc", "")
        
        # Add the grouping variable by left_join.
        df3 <- df3 %>% 
          dplyr::left_join(pData_tmp_sub %>% 
                             dplyr::select(`All ROIs`, !!as.name(grouping_var)), 
                           by = "All ROIs")
        # Replace NAs with character strings.
        df3[[grouping_var]] <- df3[[grouping_var]] %>% tidyr::replace_na('NA') # https://www.statology.org/replace-na-with-string-in-r/
        # If remove_na_imm_decon is TRUE (or true or True or tRuE or whatever), remove observations with NAs.
        if(remove_na_imm_decon | str_to_lower(remove_na_imm_decon)=="true") {
          df3 <- df3 %>% dplyr::filter(!!as.name(grouping_var) != "NA")
        }
        
        # If the number of unique values of the pData column `grouping_var` > 50, split into multiple groups for graphing. 
        # https://forum.posit.co/t/diagram-overload-split-data-into-multiple-charts/104355
        unique_values_max <- ifelse(method %in% c("quantiseq", "epic", "cibersort_abs"), 50, 5)
        unique_values <- unique(df3[[grouping_var]])
        group <- (1:length(unique_values) %/% unique_values_max) + 1
        names(group) <- unique_values
        df3$ChunkingGroup <- group[df3[[grouping_var]]]
        
        # Make sure we have enough colors.
        n_colors <- df3$cell_type %>% unique %>% length
        mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
        
        # For quantiseq, CIBERSORT absolute (when installed), and EPIC, create a stacked bar plot to show between-cell-type (within-sample) comparisons.
        if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
          # Stacked bar chart.
          # https://stackoverflow.com/questions/40361800/r-ggplot-stacked-geom-rect
          # ^ for stacked bar charts using geom_rect().
          # See also https://stackoverflow.com/questions/28956442/automatically-resize-bars-in-ggplot-for-uniformity-across-several-graphs-r
          for(group in unique(df3$ChunkingGroup)) {
            message(glue::glue("Creating stacked bar graph for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            plot <- df3 %>% 
              dplyr::filter(ChunkingGroup==group) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = fraction, fill = cell_type)) + 
              geom_bar(stat = "identity") + scale_fill_manual(values = mycolors) + # , guide = FALSE
              scale_color_manual(values = mycolors) + # , guide = FALSE
              theme_bw() + 
              theme(panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) + 
              scale_x_discrete(limits = rev(levels(imm_decon_res_list[[method]]))) + 
              labs(y = "quantity",
                title = paste0(method, " deconvolution | subset by ", subset_tag, 
                                  " | level ", subset_var_level, 
                                  "\n compartmentalized by ", grouping_var,
                                  " | group ", group))
            if(length(unique_values) > 10) {
              plot <- plot + coord_flip()
            }
            # Add to plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_stacked-bar-graph_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
            }
            
          }
        }
        
        # For CIBERSORT (when installed), MCPcounter, xCell, Abis, and Estimate, create a dot plot to show between-sample (within-cell-type) comparisons.
        if(method %in% c("cibersort", "mcp_counter", "xcell", "abis", "estimate")) {
          for(group in unique(df3$ChunkingGroup)) {
            message(glue::glue("Creating dotplot for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            plot <- df3 %>% 
              dplyr::filter(ChunkingGroup==group) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = score, color = cell_type)) +
              geom_point(size = 4) +
              facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
              scale_fill_manual(values = mycolors) + # , guide = FALSE
              scale_color_manual(values = mycolors) + # , guide = FALSE
              coord_flip() +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) +
              labs(title = paste0(method, " deconvolution | subset by ", subset_tag, 
                                  " | level ", subset_var_level, 
                                  "\n compartmentalized by ", grouping_var,
                                  " | group ", group))
            # Add to the plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_dotplot_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
            }
            
          }
        }
        
        # For SpatialDecon and MCPCounter, create a box/violin/dot plot to show between-sample (within-cell-type) comparisons
        if(method %in% c("spatialdecon", "mcp_counter")) {
          for(group in unique(df3$ChunkingGroup)) {
            message(glue::glue("Creating boxplot for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            
            # Subset to include only the cells in `spatial_decon_cells_plot`, if applicable
            cells_filter <- df3$cell_type %>% unique
            if(!flagVariable(spatial_decon_cells_plot) & (sum(spatial_decon_cells_plot %in% (df3$cell_type %>% unique)) >= 1)) cells_filter <- intersect(spatial_decon_cells_plot, df3$cell_type %>% unique)
            
            plot <- df3 %>% 
              dplyr::filter(ChunkingGroup==group & cell_type %in% cells_filter) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = score)) +
              geom_boxplot(aes(fill = !!as.name(grouping_var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
              geom_jitter(aes(color = !!as.name(grouping_var)), width = 0.15, size = 2, alpha = 0.9) + 
              labs(x = NULL, y = "Abundance score", fill = grouping_var, color = grouping_var) +
              facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
              scale_fill_manual(values = c("#FD6563","#3767A9", "grey")) + # , guide = FALSE
              scale_color_manual(values = c("#FD6563","#3767A9", "grey")) + # , guide = FALSE
              # coord_flip() +
              theme_bw() +
              theme(axis.text.x = element_blank(), # 
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) # +
              # labs(title = paste0(method, " deconvolution | subset by ", subset_tag, 
              #                     " | level ", subset_var_level, 
              #                     "\n compartmentalized by ", grouping_var,
              #                     " | group ", group))
            # Add to the plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_boxplot_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
            }
            
          }
        }
        
        # Clean up
        rm(df2)
        gc()
        
      } # End for() loop: grouping variables. 
      
    } # End for() loop: subset variable levels.
    
  } # End for() loop: subset variables.
  
} # End for() loop: deconvolution methods.

rm(pData_tmp)
gc()

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Export deconvolution results as RDS file and Microsoft Excel file. 
saveRDS(imm_decon_res_list, paste0(output_dir_rdata, "immune-deconvolution_results.rds"))
openxlsx::write.xlsx(imm_decon_res_list, file = paste0(output_dir_tabular, "immune-deconvolution_results_by-method.xlsx"))
# Export differential abundance results as RDS file and Microsoft Excel file.
if(exists("da_res_df")) {
  saveRDS(da_res_df, paste0(output_dir_rdata, "immune-deconvolution_differential-abundance_results.rds"))
  openxlsx::write.xlsx(da_res_df, file = paste0(output_dir_tabular, "immune-deconvolution_differential-abundance_results.xlsx"))
}
# Export the raw plots as RDS file.
plot_list %>% saveRDS(paste0(output_dir_rdata, "immune-deconvolution_plots-list.rds"))