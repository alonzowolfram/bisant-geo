## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
target_data_object_list <- readRDS(cl_args[5])

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

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## TCR analysis ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if(!flagVariable(module_tcr) && module_tcr %in% names(target_data_object_list)) { # Only run the module if the TCR module is provided and the combined module is in the target data object list
  # Get the TCR object
  target_data_object <- target_data_object_list[[module_tcr]]
  
  ### ................................................
  ###
  ### Analysis ----
  ###
  ### ................................................
  # Get TCR probes
  tcr_probes <- fData(target_data_object)$TargetName[base::grepl("TR[A/B/D/G][C/J/V]", fData(target_data_object)$TargetName)]
  if(length(tcr_probes) < 1) {
    warning("There are no TCR probes in the QC-ed data set. No analysis will be performed")
    plot_list <- list()
    plot_grid_list <- list()
    anova_list <- list()
  } else {
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
    
    #### ................................................
    ####
    #### Grouped analysis ----
    ####
    #### ................................................
    if(is.null(grouping_vars_tcr) || sum(grouping_vars_tcr != "") < 1) {
      pData(target_data_object)[["Full data set"]] <- "Full data set"
      grouping_vars_tcr <- c("Full data set") 
    }
    
    pdata <- pData(target_data_object)
    
    plot_list <- list()
    anova_list <- list()
    for(var in grouping_vars_tcr) {
      plot_list[[var]] <- list()
      
      # Convert the current grouping variable to factor
      pdata[[var]] <- as.factor(pdata[[var]])
      
      # Check that we have enough levels to perform ANOVA
      n_grouping_var_levels <- pdata[[var]] %>% unique %>% length
      if(n_grouping_var_levels < 2) {
        warning(glue::glue("Grouping variable {grouping_var} has fewer than 2 levels. Skipping ANOVA"))
        anova_list[[var]][["ANOVA"]][["Shannon"]] <- NA
        anova_list[[var]][["ANOVA"]][["Simpson"]] <- NA
        anova_list[[var]][["ANOVA"]][["InvSimpson"]] <- NA
        anova_list[[var]][["ANOVA"]][["Gini"]] <- NA
        anova_list[[var]][["TukeyHSD"]][["Shannon"]] <- NA
        anova_list[[var]][["TukeyHSD"]][["Simpson"]] <- NA
        anova_list[[var]][["TukeyHSD"]][["InvSimpson"]] <- NA
        anova_list[[var]][["TukeyHSD"]][["Gini"]] <- NA
      } else {
        # Perform ANOVA.
        anova_res_shannon <- aov(data = pdata, ShannonH ~ get(var), na.action=na.exclude)
        anova_res_simpson <- aov(data = pdata, Simpson ~ get(var), na.action=na.exclude)
        anova_res_invsimpson <- aov(data = pdata, InvSimpson ~ get(var), na.action=na.exclude)
        anova_res_gini <- aov(data = pdata, Gini ~ get(var), na.action=na.exclude)
        
        tukey_res_shannon <- TukeyHSD(anova_res_shannon)
        tukey_res_simpson <- TukeyHSD(anova_res_simpson)
        tukey_res_invsimpson <- TukeyHSD(anova_res_invsimpson)
        tukey_res_gini <- TukeyHSD(anova_res_gini)
        
        anova_list[[var]][["ANOVA"]][["Shannon"]] <- anova_res_shannon
        anova_list[[var]][["ANOVA"]][["Simpson"]] <- anova_res_simpson
        anova_list[[var]][["ANOVA"]][["InvSimpson"]] <- anova_res_invsimpson
        anova_list[[var]][["ANOVA"]][["Gini"]] <- anova_res_gini
        anova_list[[var]][["TukeyHSD"]][["Shannon"]] <- tukey_res_shannon
        anova_list[[var]][["TukeyHSD"]][["Simpson"]] <- tukey_res_simpson
        anova_list[[var]][["TukeyHSD"]][["InvSimpson"]] <- tukey_res_invsimpson
        anova_list[[var]][["TukeyHSD"]][["Gini"]] <- tukey_res_gini
      }
      
      # Make sure we have enough colors
      n_colors <- pdata[[var]] %>% unique %>% length
      mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
      
      # Graph
      plot_shannon <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = ShannonH, fill = !!as.name(var))) + 
        geom_boxplot(aes(fill = !!as.name(var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
        geom_jitter(aes(color = !!as.name(var)), width = 0.15, size = 2, alpha = 0.9) + 
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        theme_bw() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
        ) +
        labs(title = glue::glue(""),
             x = paste0(""),
             y = paste0("Shannon"))
      plot_simpson <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = Simpson, fill = !!as.name(var))) + 
        geom_boxplot(aes(fill = !!as.name(var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
        geom_jitter(aes(color = !!as.name(var)), width = 0.15, size = 2, alpha = 0.9) + 
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        theme_bw() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
        ) +
        labs(title = glue::glue(""),
             x = paste0(""),
             y = paste0("Simpson"))
      plot_invsimpson <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = InvSimpson, fill = !!as.name(var))) + 
        geom_boxplot(aes(fill = !!as.name(var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
        geom_jitter(aes(color = !!as.name(var)), width = 0.15, size = 2, alpha = 0.9) + 
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        theme_bw() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
        ) +
        labs(title = glue::glue(""),
             x = paste0(""),
             y = paste0("InvSimpson"))
      plot_gini <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = Gini, fill = !!as.name(var))) + 
        geom_boxplot(aes(fill = !!as.name(var)), width = 0.6, alpha = 0.7, outlier.shape = NA) + 
        geom_jitter(aes(color = !!as.name(var)), width = 0.15, size = 2, alpha = 0.9) + 
        scale_fill_manual(values = mycolors) + # , guide = FALSE
        scale_color_manual(values = mycolors) + # , guide = FALSE
        theme_bw() + 
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right",
          strip.background = element_rect(fill = "#E4E4E4", color = "#333333")
        ) +
        labs(title = glue::glue(""),
             x = paste0(""),
             y = paste0("Gini"))
      # Add to plot list
      plot_list[[var]][["Shannon"]] <- plot_shannon
      plot_list[[var]][["Simpson"]] <- plot_simpson
      plot_list[[var]][["InvSimpson"]] <- plot_invsimpson
      plot_list[[var]][["Gini"]] <- plot_gini
    }
    
    # Arrange graphs
    plot_grid_list <- list()
    for(var in names(plot_list)) {
      # Get the list of plots for variable `var`
      p_list <- plot_list[[var]]
      
      # Get the number of levels in the grouping variable
      n_levels <- pdata[[var]] %>% unique %>% length
      # We will use the number of levels to determine whether we arrange the 3 diversity graphs horizontally or vertically
      nCol <- 2 #ifelse(n_levels < 4, 3, 1)
      
      # Extract legend from one plot
      legend <- ggpubr::get_legend(p_list[[1]] + 
                                      guides(color = "none") + 
                                      theme(legend.position = "right", legend.box.margin = margin(0,0,0,0)))
      # Strip legends from p_list
      for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}
      
      # Arrange plots in p_list onto a grid
      plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
      plot_grid <- plot_grid %>% ggpubr::as_ggplot() %>% 
        ggpubr::annotate_figure(left = grid::textGrob("", 
                                                      hjust = 0, 
                                                      rot = 90, 
                                                      vjust = 1#, 
                                                      # gp = grid::gpar(cex = scaling_factor)
        ),
        bottom = grid::textGrob(""#, 
                                #gp = grid::gpar(cex = scaling_factor)
        ),
        top = grid::textGrob(glue::glue("TCR diversity and distribution | grouping variable: {var}")#, 
                             # gp = grid::gpar(cex = scaling_factor)
        )) 
      # Add the legend
      final_grob <- gridExtra::arrangeGrob(plot_grid, legend, ncol = 2, widths = c(1, 0.25))
      
      # Add to the plot list
      plot_grid_list[[var]] <- final_grob
      # Save plot to disk
      final_grob
      file_name <- glue::glue("TCR_plots_by_{var}")
      for(file_type in output_plot_file_types) {
        ggsave(
          plot = final_grob,
          filename = glue::glue("{file_name}.{file_type}"), 
          path = output_dir_imgs, 
          width = 12, height = 9, units = "in")
      }
      
    } # End for() loop: grouping variables
    
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

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export NanoStringGeoMxSet as RDS file
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_TCR-analysis.rds"))
# Export the raw plots as RDS file
if(exists("plot_list")) plot_list %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_plots-list.rds"))
if(exists("plot_grid_list")) plot_grid_list %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_plot-grid-list.rds"))
# Export ANOVA results
if(exists("anova_list")) saveRDS(anova_list, paste0(output_dir_rdata, "TCR-analysis_ANOVA-res-list.rds"))
