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

###################################################################
##
## TCR analysis
##
###################################################################
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

if(!(is.null(module_tcr) || module_tcr == "")) { # Only run the module if the TCR module is provided.
  # Add a section header.
  pptx <- pptx %>%
    officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
    officer::ph_with(value = paste0("TCR analysis"),
                     location = ph_location_label(ph_label = "Title 1"))
  
  ## ----------------------------------------------------------------
  ##
  ## Subsetting
  ##
  ## ----------------------------------------------------------------
  # Subset TCR + WTA prior to normalization.
  target_data_object_tcr <- subset(target_data_object, Module %in% module_tcr)
  
  ## ----------------------------------------------------------------
  ##
  ## Normalization
  ##
  ## ----------------------------------------------------------------
  # Subtract background
  target_data_object_tcr <- normalize(target_data_object_tcr,
                                          norm_method = "subtractBackground",
                                          fromElt = "exprs",
                                          toElt = "bgsub"
  )
  
  # Q3 normalization
  target_data_object_tcr <- normalize(target_data_object_tcr,
                                          norm_method = "quant",
                                          desiredQuantile = 0.9,
                                          fromElt = "bgsub",
                                          toElt = "q_norm_bgsub"
  )
  
  ## ----------------------------------------------------------------
  ##
  ## Analysis
  ##
  ## ----------------------------------------------------------------
  # Get TCR probes
  tcr_probes <- fData(target_data_object_tcr)$TargetName[base::grepl("TR[A/B/D/G][C/J/V]", fData(target_data_object_tcr)$TargetName)]
  if(length(tcr_probes) < 1) {
    warning("There are no TCR probes in the QC-ed data set. No analysis will be performed.")
    plot_list <- list()
    plot_grid_list <- list()
    anova_list <- list()
  } else {
    # Calculate distribution (Gini coefficient).
    gini <- apply(assayDataElement(target_data_object_tcr, elt = "bgsub")[tcr_probes, ], 2, calculate_gini_coefficient)
    gini[!is.finite(gini)] <- NA
    pData(target_data_object)$Gini <- gini
    
    # Calculate diversity (Shannon, Simpson, and inverse Simpson).
    shannon_h <- vegan::diversity(t(assayData(target_data_object_tcr)$bgsub[tcr_probes, ]), index = "shannon", MARGIN = 1)
    simpson <- vegan::diversity(t(assayData(target_data_object_tcr)$bgsub[tcr_probes, ]), index = "simpson", MARGIN = 1)
    invsimpson <- vegan::diversity(t(assayData(target_data_object_tcr)$bgsub[tcr_probes, ]), index = "invsimpson", MARGIN = 1)
    shannon_h[!is.finite(shannon_h)] <- NA
    simpson[!is.finite(simpson)] <- NA
    invsimpson[!is.finite(invsimpson)] <- NA
    pData(target_data_object)$ShannonH <- shannon_h
    pData(target_data_object)$Simpson <- simpson
    pData(target_data_object)$InvSimpson <- invsimpson
    
    ## ----------------------------------------------------------------
    ##
    ## Grouped analysis
    ##
    ## ----------------------------------------------------------------
    if(is.null(tcr_grouping_vars) || sum(tcr_grouping_vars != "") < 1) {
      pData(target_data_object)[["Full data set"]] <- "Full data set"
      tcr_grouping_vars <- c("Full data set") 
    }
    
    pdata <- pData(target_data_object)
    
    plot_list <- list()
    anova_list <- list()
    for(var in tcr_grouping_vars) {
      plot_list[[var]] <- list()
      
      # Convert the current grouping variable to factor.
      pdata[[var]] <- as.factor(pdata[[var]])
      
      # Check that we have enough levels to perform ANOVA.
      n_grouping_var_levels <- pdata[[var]] %>% unique %>% length
      if(n_grouping_var_levels < 2) {
        warning("Grouping variable ", grouping_var, " has fewer than 2 levels. Skipping ANOVA.")
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
      
      # Make sure we have enough colors.
      n_colors <- pdata[[var]] %>% unique %>% length
      mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
      
      # Graph.
      plot_shannon <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = ShannonH, color = !!as.name(var), fill = !!as.name(var))) + 
        geom_boxplot() + 
        theme_bw()
      plot_simpson <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = Simpson, color = !!as.name(var), fill = !!as.name(var))) + 
        geom_boxplot() + 
        theme_bw()
      plot_invsimpson <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = InvSimpson, color = !!as.name(var), fill = !!as.name(var))) + 
        geom_boxplot() + 
        theme_bw()
      plot_gini <- pdata %>% 
        ggplot(aes(x = !!as.name(var), y = Gini, color = !!as.name(var), fill = !!as.name(var))) + 
        geom_boxplot() + 
        theme_bw()
      # Add to plot list.
      plot_list[[var]][["Shannon"]] <- plot_shannon
      plot_list[[var]][["Simpson"]] <- plot_simpson
      plot_list[[var]][["InvSimpson"]] <- plot_invsimpson
      plot_list[[var]][["Gini"]] <- plot_gini
    }
    
    # Arrange graphs and add to PowerPoint.
    # Graphing parameters.
    plot_width <- 12
    plot_height <- 12 
    units <- "in"
    res <- 300
    plot_grid_list <- list()
    for(var in names(plot_list)) {
      # Get the list of plots for variable `var`.
      p_list <- plot_list[[var]]
      
      # Get the number of levels in the grouping variable. 
      n_levels <- pdata[[var]] %>% unique %>% length
      # We will use the number of levels to determine whether we arrange the 3 diversity graphs horizontally or vertically.
      nCol <- 2 #ifelse(n_levels < 4, 3, 1)
      
      # Strip legends from p_list.
      for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}
      
      # Arrange plots in p_list onto a grid.
      plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
      plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("", 
                                                                               hjust = 0, 
                                                                               rot = 90, 
                                                                               vjust = 1#, 
                                                                               # gp = grid::gpar(cex = scaling_factor)
      ),
      bottom = grid::textGrob(""#, 
                              #gp = grid::gpar(cex = scaling_factor)
      ),
      top = grid::textGrob(paste0("TCR diversity and distribution | grouping variable: ", var)#, 
                           # gp = grid::gpar(cex = scaling_factor)
      ))
      
      # Save to list.
      plot_grid_list[[var]] <- plot_grid
      
      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, 
                         paste0("TCR-diversity-distribution_graphs_by-", var, ".eps") %>% regexPipes::gsub("\\/", "_"))
      png_path <- paste0(output_dir_pubs, 
                         paste0("TCR-diversity-distribution_graphs_by-", var, ".png") %>% regexPipes::gsub("\\/", "_"))
      saveEPS(plot_grid, eps_path, width = (plot_width), height = (plot_height))
      savePNG(plot_grid, png_path, width = (plot_width), height = (plot_height), units = units, res = (res))
      
      # Add to the PowerPoint.
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("TCR diversity and distribution"),
                         location = ph_location_label(ph_label = "Title 1")) %>%
        officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                         location = ph_location_label(ph_label = "Content Placeholder 2"),
                         use_loc_size = FALSE)
    } # End for() loop: grouping variables.
  } # End if() there are TCR probes.

} # End if() the TCR module is provided.

###################################################################
##
## Export to disk
##
###################################################################
# Export PowerPoint file.
print(pptx, cl_args[5])
# Export NanoStringGeoMxSet as RDS file.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_TCR-analysis.rds"))
# Export the raw plots as RDS file.
if(exists("plot_list")) plot_list %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_plots-list.rds"))
if(exists("plot_grid_list")) plot_grid_list %>% saveRDS(paste0(output_dir_rdata, "TCR-analysis_plot-grid-list.rds"))
# Export ANOVA results.
if(exists("anova_list")) saveRDS(anova_list, paste0(output_dir_rdata, "TCR-analysis_ANOVA-res-list.rds"))
