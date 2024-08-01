###################################################################
##                                                                
## Setup 
##
###################################################################
## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object.
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

# Set the normalization method.
normalization_method <- normalization_names[names(normalization_names)==normalization_methods[1]]

# Cell profile matrices for SpatialDecon.
data("safeTME")
data("safeTME.matches")

# # Set path to CIBERSORT required files.
# set_cibersort_binary(path_to_cibersort)
# set_cibersort_mat(path_to_lm22)
# 
# # Calculate TPM - this is necessary for CIBERSORT among others, not so much for xCell or MCP-counter.
# # https://bioinformatics.stackexchange.com/questions/2567/how-can-i-calculate-gene-length-for-rpkm-calculation-from-counts-data
# # But can we even calculate TPM for GeoMx data? Bc it's probe-based, so it wouldn't have the same assumptions that RNA-seq does ... 

# Add a section header.
pptx <- pptx %>%
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Immune deconvolution"),
                   location = ph_location_label(ph_label = "Title 1"))

###################################################################
##
## Immune deconvolution
##
###################################################################
## ----------------------------------------------------------------
##
## Calculation
##
## ----------------------------------------------------------------
# Subset the target data object to include only the transcriptome probes.
target_data_object_exprs <- subset(target_data_object, Module %in% modules_exprs)
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
    warning(paste0(method, " is currently not supported by this pipeline. Please note that TIMER and ConsensusTME, while included in the immunedeconv package, are currently not available in this pipeline due to extra arguments that must be passed to the function; and CIBERSORT will not be available until we figure out how to make the source code play nicely with the immunedeconv package."))
    next
  } else {
    # Error catching: https://stackoverflow.com/a/55937737/23532435
    skip_to_next <- FALSE
    
    # If the method is mmcp_counter, check if the species is mouse (Mus musculus).
    if(method == "mmcp_counter") {
      if(species != "Mus musculus") {
        warning(paste0("mMCP-counter can only be used with human data. Skipping to the next one."))
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
    
    # Special steps needed for SpatialDecon.
    # Also, we can only run it on human data?
    if(method == "spatialdecon") {
      # The spatialdecon function takes 3 arguments of expression data:
      #   
      # 1. The normalized data.
      # 2. A matrix of expected background for all data points in the normalized data matrix.
      # 3. Optionally, either a matrix of per-data-point weights, or the raw data, which is used to derive weights (low counts are less statistically stable, and this allows spatialdecon to down-weight them.)
      # https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html
      
      # Estimate each data point's expected BG from the negative control probes from its corresponding observation.
      negnames <- fData(target_data_object_exprs) %>% 
        dplyr::filter(Negative == TRUE & Module %in% modules_exprs) %>% 
        .$TargetName
      bg <- derive_GeoMx_background(norm = target_data_object_exprs@assayData[[normalization_methods[1]]],
                                    probepool = fData(target_data_object_exprs)$Module,
                                    negnames = negnames)
      
      signif(safeTME[seq_len(3), seq_len(3)], 2)
      # heatmap(sweep(safeTME, 1, apply(safeTME, 1, max), "/"),
      #         labRow = NA, margins = c(10, 5))
      
      # Set the cell profile matrix based on the species.
      if(species == "Homo sapiens") {
        cpm <- safeTME
      } else {
        cpm <- profile_matrix %>% as.matrix
      }
      
      # Run spatial deconvolution.
      system.time({
        res <- tryCatch(runspatialdecon(object = target_data_object_exprs,
                                        norm_elt = normalization_methods[1], # "neg_norm"   "log_norm"   "bg_sub"     "exprs"      "bg_sub_neg" "quant"      "bg_sub_q3"  "q3_norm"   
                                        raw_elt = "exprs",
                                        X = cpm,
                                        align_genes = TRUE),
                        error = function(e) {skip_to_next <<- TRUE})
      })
      if(class(res) != "logical") imm_decon_res <- res$beta %>% t %>% as.data.frame %>% rownames_to_column("cell_type")
    } else {
      imm_decon_res <- tryCatch(immunedeconv::deconvolute(exprs_mat_effective, method),
                                error = function(e) {skip_to_next <<- TRUE})
    }
    
    if(skip_to_next) {
      warning(paste0("An error occurred when trying to run immune deconvolution method ", method, ". Skipping to the next method."))
      next
    }
    
    imm_decon_res_list[[method]] <- imm_decon_res
  }
}

## ----------------------------------------------------------------
##
## Visualization
##
## ----------------------------------------------------------------
# https://omnideconv.org/immunedeconv/articles/detailed_example.html
# Parameters.
plot_width <- 12
plot_height <- 12 
units <- "in"
res <- 300
# Rename the samples to reflect their segment and type.
pData_tmp <- pData(target_data_object) %>% as.data.frame %>% tibble::rownames_to_column(var = "rownames")
observation_identifiers <- intersect(observation_identifiers, colnames(pData_tmp)) # Make sure the observation identifiers are actually in the pData.
pData_tmp <- pData_tmp %>% tidyr::unite("All ROIs", c(observation_identifiers, rownames), remove = FALSE, na.rm = FALSE, sep = " | ") 
pData_tmp$`All ROIs` <- pData_tmp$`All ROIs` %>% regexPipes::gsub("\\.dcc", "")
pData_tmp$`Complete dataset` <- "Complete dataset"
plot_list <- list()
for(method in names(imm_decon_res_list)) {
  plot_list[[method]] <- list()
  df <- imm_decon_res_list[[method]]

  # QuanTIseq, CIBERSORT (absolute), Epic, SpatialDecon - visualize as stacked bar charts.
  # MCP-counter - visualize as dot plot.
  
  for(subset_var in imm_decon_subset_vars) {
    if(is.na(subset_var) || subset_var == "NA") {
      subset_tag <- "All ROIs"
      subset_var <- "Complete dataset"
    } else {
      subset_tag <- subset_var
    }
    
    # Get the levels of the current subset_var.
    subset_var_levels <- pData_tmp[[subset_var]] %>% as.factor %>% levels # as.factor needed because it might be a character vector.
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      # Get all the samples belonging to the current subset_var_level.
      ind <- pData_tmp[[subset_var]] == subset_var_level
      pData_tmp_sub <- pData_tmp[ind,]
      
      for(grouping_var in imm_decon_grouping_vars) {
        df2 <- df
        
        # Gather the dataframe for graphing.
        if(method %in% c("quantiseq", "epic", "cibersort_abs", "spatialdecon")) {
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
        # If imm_decon_remove_na is TRUE (or true or True or tRuE or whatever), remove observations with NAs.
        if(imm_decon_remove_na | str_to_lower(imm_decon_remove_na)=="true") {
          df3 <- df3 %>% dplyr::filter(!!as.name(grouping_var) != "NA")
        }
        
        # If the number of unique values of the pData column `grouping_var` > 50, split into multiple groups for graphing. 
        # https://forum.posit.co/t/diagram-overload-split-data-into-multiple-charts/104355
        unique_values_max <- ifelse(method %in% c("quantiseq", "epic", "cibersort_abs"), 50, 5)
        unique_values <- unique(df3[[grouping_var]])
        group <- (1:length(unique_values) %/% unique_values_max) + 1
        names(group) <- unique_values
        df3$Group <- group[df3[[grouping_var]]]
        
        # Make sure we have enough colors.
        n_colors <- df3$cell_type %>% unique %>% length
        mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
        
        # For quantiseq, CIBERSORT absolute (when installed), and EPIC, create a stacked bar plot to show between-cell-type (within-sample) comparisons.
        if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
          # Stacked bar chart.
          # https://stackoverflow.com/questions/40361800/r-ggplot-stacked-geom-rect
          # ^ for stacked bar charts using geom_rect().
          # See also https://stackoverflow.com/questions/28956442/automatically-resize-bars-in-ggplot-for-uniformity-across-several-graphs-r
          for(group in unique(df3$Group)) {
            plot <- df3 %>% 
              dplyr::filter(Group==group) %>% 
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
            plot_list[[method]][[group]] <- plot
            
            # Save to EPS and PNG and then ...
            eps_path <- paste0(output_dir_pubs, 
                               paste0(method, "_subset-by-", subset_tag, "_subset-", subset_var_level, "_compartment-var-", grouping_var, "_group-", group, "_between-cell-type-comparison.eps") %>% regexPipes::gsub("\\/", "_"))
            png_path <- paste0(output_dir_pubs, 
                               paste0(method, "_subset-by-", subset_tag, "_subset-", subset_var_level, "_compartment-var-", grouping_var, "_group-", group, "_between-cell-type-comparison.png") %>% regexPipes::gsub("\\/", "_"))
            saveEPS(plot, eps_path, width = plot_width, height = plot_height)
            savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
            
            # Save to PPT if there won't be too many graphs.
            if(length(unique(df3$Group)) <= 7) {
              pptx <- pptx %>%
                officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
                officer::ph_with(value = paste0("Immune deconvolution"),
                                 location = ph_location_label(ph_label = "Title 1")) %>% 
                officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                                 location = ph_location_label(ph_label = "Content Placeholder 2"),
                                 use_loc_size = FALSE) 
            } else {
              warning(paste0("The method ", method, " will generate an overwhelming number of graphs to place in a PowerPoint, so we're going to save the graphs to disk but not into the generated PowerPoint."))
            }
          }
          
        }
        
        # For CIBERSORT (when installed), MCPcounter, xCell, Abis, and Estimate, create a dot plot to show between-sample (within-cell-type) comparisons.
        if(method %in% c("cibersort", "mcp_counter", "xcell", "abis", "estimate")) {
          for(group in unique(df3$Group)) {
            plot <- df3 %>% 
              dplyr::filter(Group==group) %>% 
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
            
            plot_list[[method]][[group]] <- plot
            
            # Save to EPS and PNG and then ...
            eps_path <- paste0(output_dir_pubs, 
                               paste0(method, "_subset-by-", subset_tag, "_subset-", subset_var_level, "_compartment-var-", grouping_var, "_group-", group, "_between-sample-comparison.eps") %>% regexPipes::gsub("\\/", "_"))
            png_path <- paste0(output_dir_pubs, 
                               paste0(method, "_subset-by-", subset_tag, "_subset-", subset_var_level, "_compartment-var-", grouping_var, "_group-", group, "_between-sample-comparison.png") %>% regexPipes::gsub("\\/", "_"))
            saveEPS(plot, eps_path, width = plot_width, height = plot_height)
            savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
            
            # Save to PPT if there won't be too many graphs.
            if(length(unique(df3$Group)) <= 5) {
              pptx <- pptx %>%
                officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
                officer::ph_with(value = paste0("Immune deconvolution"),
                                 location = ph_location_label(ph_label = "Title 1")) %>% 
                officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                                 location = ph_location_label(ph_label = "Content Placeholder 2"),
                                 use_loc_size = FALSE) 
            } else {
              warning(paste0("The method ", method, " will generate an overwhelming number of graphs to place in a PowerPoint, so we're going to save the graphs to disk but not into the generated PowerPoint."))
            }
          }
        }
        rm(df2)
        gc()
        
      } # End for() loop: grouping variables. 
      
    } # End for() loop: subset variable levels.
    
  } # End for() loop: subset variables.
  
} # End for() loop: deconvolution methods.

rm(pData_tmp)
gc()

###################################################################
##
## Export to disk
##
###################################################################

# Export PowerPoint file.
print(pptx, cl_args[5])
# Export NanoStringGeoMxSet as RDS file.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_immune-deconvolution.rds"))
# Export deconvolution results as RDS file. 
saveRDS(imm_decon_res_list, paste0(output_dir_rdata, "immune-deconvolution_results.rds"))
# Export deconvolution results as Microsoft Excel file. 
openxlsx::write.xlsx(imm_decon_res_list, file = paste0(output_dir_tabular, "immune-deconvolution_results_by-method.xlsx"))
# Export the raw plots as RDS file.
plot_list %>% saveRDS(paste0(output_dir_rdata, "immune-deconvolution_plots-list.rds"))