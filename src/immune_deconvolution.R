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
# Extract the expression matrix.
exprs_mat <- target_data_object@assayData[[normalization_methods[1]]]
# Loop through imm_decon_methods.
imm_decon_res_list <- list()
for(method in imm_decon_methods) {
  if(!(method %in% c("quantiseq", "mcp_counter", "xcell", "epic", "abis", "estimate"))) {
    warning(paste0(method, " is currently not supported by this pipeline. Please note that TIMER and ConsensusTME, while included in the immunedeconv package, are currently not available in this pipeline due to extra arguments that must be passed to the function; and CIBERSORT will not be available until we figure out how to make the source code play nicely with the immunedeconv package."))
    next
  } else {
    # Error catching: https://stackoverflow.com/a/55937737/23532435
    skip_to_next <- FALSE
    
    imm_decon_res <- tryCatch(immunedeconv::deconvolute(exprs_mat, method),
                                             error = function(e) {skip_to_next <<- TRUE})
    
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
pData_tmp <- pData(target_data_object) %>% as.data.frame
pData_tmp <- pData_tmp %>% dplyr::mutate(NovoLabel = paste0(!!as.name(individual_identifier), " | ", !!as.name(compartment_identifier), " | ", !!as.name(segment_identifier), " | ", rownames(.) %>% regexPipes::gsub("\\.dcc", "")))
plot_list <- list()
for(method in names(imm_decon_res_list)) {
  df <- imm_decon_res_list[[method]]

  # QuanTIseq, CIBERSORT (absolute) - visualize as stacked bar charts.
  # MCP-counter - visualize as ??
  
  # Rename the samples to reflect their segment and type.
  if(identical(colnames(df)[2:ncol(df)], rownames(pData_tmp))) colnames(df)[2:ncol(df)] <- pData_tmp$NovoLabel
  
  # Gather the dataframe for graphing.
  if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
    df2 <- df %>% 
      gather(sample, fraction, -cell_type)
  } else {
    df2 <- df %>% 
      gather(sample, score, -cell_type)
  }
  
  # If the number of samples > 50, split into multiple groups for graphing. 
  # https://forum.posit.co/t/diagram-overload-split-data-into-multiple-charts/104355
  sample_max <- ifelse(method %in% c("quantiseq", "epic", "cibersort_abs"), 50, 5)
  samples <- unique(df2$sample)
  group <- (1:length(samples) %/% sample_max) + 1
  names(group) <- samples
  df2$Group <- group[df2$sample]
  
  # Make sure we have enough colors.
  n_colors <- df2$cell_type %>% unique %>% length
  mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
  
  # For quantiseq, CIBERSORT absolute (when installed), and EPIC, create a stacked bar plot to show between-cell-type (within-sample) comparisons.
  if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
    # Stacked bar chart.
    # https://stackoverflow.com/questions/40361800/r-ggplot-stacked-geom-rect
    # ^ for stacked bar charts using geom_rect().
    # See also https://stackoverflow.com/questions/28956442/automatically-resize-bars-in-ggplot-for-uniformity-across-several-graphs-r
    for(group in unique(df2$Group)) {
      plot <- df2 %>% 
        dplyr::filter(Group==group) %>% 
        ggplot(aes(x = sample, y = fraction, fill = cell_type)) + 
        geom_bar(stat = "identity") + 
        coord_flip() + 
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        theme_bw() + 
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) + 
        scale_x_discrete(limits = rev(levels(imm_decon_res_list[[method]]))) + 
        labs(title = paste0(method, " deconvolution | group ", group))
      plot_list[[method]][[group]] <- plot
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, method, "_group-", group, "_between-cell-type-comparison.eps")
      png_path <- paste0(output_dir_pubs, method, "_group-", group, "_between-cell-type-comparison.png")
      saveEPS(plot, eps_path, width = plot_width, height = plot_height)
      savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
      
      # Save to PPT if there won't be too many graphs.
      if(length(unique(df2$Group)) <= 7) {
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
  
  # For quantiseq, CIBERSORT absolute (when installed), and EPIC, create a dot plot to show between-sample (within-cell-type) comparisons.
  if(method %in% c("cibersort", "mcp_counter", "xcell", "abis", "estimate")) {
    for(group in unique(df2$Group)) {
      plot <- df2 %>% 
        dplyr::filter(Group==group) %>% 
        ggplot(aes(x = sample, y = score, color = cell_type)) +
        geom_point(size = 4) +
        facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        labs(title = paste0(method, " deconvolution | group ", group))
          
      plot_list[[method]][[group]] <- plot
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, method, "_group-", group, "_between-sample-comparison.eps")
      png_path <- paste0(output_dir_pubs, method, "_group-", group, "_between-sample-comparison.png")
      saveEPS(plot, eps_path, width = plot_width, height = plot_height)
      savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
      
      # Save to PPT if there won't be too many graphs.
      if(length(unique(df2$Group)) <= 5) {
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
  rm(df)
  gc()
}

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