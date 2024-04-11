## Source the setup.R file.
source("src/setup.R")

## ---------------------------
## Test values
# pathway_table_file <- NULL
# species <- "Homo sapiens"
# individual_pathways <- c("HALLMARK_ADIPOGENESIS", "REACTOME_SHC1_EVENTS_IN_ERBB2_SIGNALING", "REACTOME_SEROTONIN_RECEPTORS", "BIOCARTA_PEPI_PATHWAY")
# subset_var <- unique(results2$`Subset variable`)[1]
# subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
# subset_var_level <- subset_var_levels[1]
# model_numbers <- results2$`Model number` %>% unique
# model_number <- model_numbers[1]
# contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model number`==model_number) %>% .$Contrast %>% unique
# contrast <- contrasts[1]
# results2 <- read.csv("/Users/lwfong/Library/CloudStorage/OneDrive-InsideMDAnderson/PRIME-TR/LFI-001_bisantine-geo/out/DSP-018-SBNET_differential-expression_2024-4-10-20-49-19/tabular/LMM-differential-expression_results.csv", row.names = 1, check.names = FALSE)

## ---------------------------
# Setup

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])
# Read in the DE genes table.
results2 <- read.csv(cl_args[6], row.names = 1, check.names = FALSE)
print(paste0("results2 dimensions: ", dim(results2)))

# If pathway_table_file is provided, check that it's a valid file.
# If not, default to hallmark, BioCarta, and Reactome pathways.
# In future iterations, we'll probably include a default pathway_table_file instead of hard-coding the pathways here.
set.default.pathways <- function() {
  cats <- c("H", rep("C2", 2))
  subcats <- c(NA, "CP:BIOCARTA", "CP:REACTOME")
  return(list(cats = cats, subcats = subcats))
}

if(!is.null(pathway_table_file)) { 
  if(!file.exists(pathway_table_file)) {
    msigdb_list <- set.default.pathways()
  } else {
    # Read in the file and check that it has at least one entry. 
    if(base::grepl("\\.csv$", pathway_table_file)) { pathway_table <- read.csv(pathway_table_file, header = F)}
    else if(base::grepl("\\.tsv$", pathway_table_file)) { pathway_table <- read.table(pathway_table_file, header = F, sep = "\t") }
    else if(base::grepl("\\.xls.*$", pathway_table_file)) { pathway_table <- read_excel(pathway_table_file, col_names = F)} 
    else {
      warning("The pathway table provided is not in CSV, TSV, or Excel format. Check the file extension. Using default pathways instead.")
      pathway_table <- NULL
      msigdb_list <- set.default.pathways()}
  }
} else {
  msigdb_list <- set.default.pathways()
}
if(exists("pathway_table")) if(!is.null(pathway_table)) msigdb_list <- as.list(pathway_table); names(msigdb_list) <- c("cats", "subcats")
  
# Check that species is valid.
if(!(species %in% c("Homo sapiens", "Mus musculus"))) stop("Please provide a valid species: either 'Homo sapiens' or 'Mus musculus.'")

## ---------------------------
# Add section header to PPT.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Pathway analysis (FGSEA)"), 
                   location = ph_location_label(ph_label = "Title 1"))

## ---------------------------
# Run pathway analysis (FGSEA).

# Load the gene sets (from MSigDB).
# https://rpubs.com/LiYumei/806213
gene_sets <- data.frame()
for(i in 1:length(msigdb_list[["subcats"]])) {
  subcat <- msigdb_list[["subcats"]][i]
  cat <- msigdb_list[["cats"]][i]
  if(is.na(subcat)) {
    subcat <- NULL
  } else {
    subcat <- subcat
  }
  
  gene_sets <- rbind(gene_sets, 
                     msigdbr(species = species, category = cat, subcategory = subcat) %>% 
                       dplyr::distinct(gs_name, gene_symbol) %>% 
                       as.data.frame()
  )
}

# Subset to include only the pathways of interest if individual_pathways is provided.
if(sum(!is.null(individual_pathways)) > 0 & sum(individual_pathways != "") > 0) gene_sets <- gene_sets %>% dplyr::filter(gs_name %in% individual_pathways)

# Create an MSigDB pathway list.
msigdbr_pathway_list <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

# Loop through DE gene table and perform FGSEA.
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
plot_list_pathway_analysis <- list()
pathway_df <- data.frame()
model_numbers <- results2$`Model number` %>% unique
for(subset_var in unique(results2$`Subset variable`)) { # We're not naming it subset_vars because we already have a variable by that name. ... 
  
  subset_var_levels <- results2 %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model_number in model_numbers) {
      i <- model_number %>% regexPipes::gsub("model_", "")
      
      contrasts <- results2 %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model number`==model_number) %>% .$Contrast %>% unique
      for(contrast in contrasts) {
        # Subset. 
        results2_sub <- results2 %>% dplyr::filter(
          `Subset variable`==subset_var &
            `Subset level`==subset_var_level &
            `Model number`==i &
            Contrast==contrast
        )
        
        # Extract the DEG data--these will be the ranks we use as input to FGSEA. 
        ranks <- results2_sub$Estimate
        names(ranks) <- results2_sub$Gene
        ranks <- ranks %>% .[order(.)]
        
        # Run FGSEA.
        fgseaRes <- fgsea(pathways = msigdbr_pathway_list, 
                          stats    = ranks,
                          minSize  = 1, # 15
                          maxSize  = Inf) # 500
        df_sub <- fgseaRes %>% dplyr::filter(is.finite(NES))
        
        # Add pathway ranking score (-log10(padj) * |NES|)
        df_sub$PathwayScore <- -log10(df_sub$padj) * abs(df_sub$NES)
        # Before calculating percentiles for pathway, cull pathway list down to its final form (i.e., what will actually be graphed.)
        # If individual_pathways is set, subset to include only the pathways of interest.
        if(!is.null(individual_pathways) & sum(individual_pathways=="") < length(individual_pathways)) df_sub <- df_sub %>% dplyr::filter(pathway %in% individual_pathways)
        # If the number of pathways > limits for graphing, perform cutoff by PathwayScore.
        if(!is.null(n_max_pathways) & n_max_pathways != "") if(nrow(df_sub) > n_max_pathways) df_sub <- df_sub %>% dplyr::top_n(n = n_max_pathways, wt = PathwayScore)
        
        # Add column indicating percentile (how high up the list a given pathway is.)
        df_sub <- df_sub %>% dplyr::arrange(-PathwayScore) %>% dplyr::mutate(percentile = 100 - (1:nrow(.) / nrow(.) * 100))
        df_sub$percentile <- format(round(df_sub$percentile, 1), nsmall = 1)
        
        # Graph bar charts.
        upper_limit <- ifelse(max(df_sub$NES) > 0, max(df_sub$NES) %>% ceiling(), max(df_sub$NES) %>% abs %>% ceiling())
        lower_limit <- ifelse(min(df_sub$NES) < 0, min(df_sub$NES) %>% floor(), min(df_sub$NES) %>% -. %>% floor())
        upper_limit_axis <- ifelse(max(df_sub$NES) < 0, 0, upper_limit)
        lower_limit_axis <- ifelse(min(df_sub$NES) > 0, 0, lower_limit)
        nes_palette <- colorRampPalette(c("blue", "white", "red"))(100)
        dodgewidth <- position_dodge(width=0.9)
        
        # Get the contrast elements for the plot label.
        contrast_element_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1]
        contrast_element_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2]
        # Graph.
        plot <- ggplot(df_sub, aes(x = NES, y = reorder(pathway, NES), fill = NES)) +
          # Bar graph.
          geom_bar(stat = "identity", color = "black") +
          # Set the colors.
          scale_fill_gradientn(colors = nes_palette, limits = c(lower_limit, upper_limit)) +
          # Change the axis limits.
          xlim(lower_limit_axis, upper_limit_axis) + 
          # # Add annotation (percentiles).
          # geom_text(aes(label = percentile), hjust = -0.25) + 
          # Add labels.
          labs(
              x = paste0(contrast_element_2, " <-> ", contrast_element_1), 
               # x = "Normalized enrichment score", 
               y = ""
               ) +
          # Use B+W theme.
          theme_bw() +
          # Position the legend on the right.
          theme(legend.position = "right")
        
        # Add to the list.
        plot_list_pathway_analysis[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]] <- plot
        
        # Add df_sub to pathway_df
        pathway_df <- rbind(pathway_df, df_sub %>% dplyr::select(-leadingEdge))

        # Graph enrichment plots.
        for(pathway in names(msigdbr_pathway_list)) {
          plotEnrichment(msigdbr_pathway_list[[pathway]], ranks, gseaParam = 1, ticksSize = 0.2) +
            ggtitle(paste0(pathway, " | contrast: ", contrast, " | subset by ", subset_var, " | level: ", subset_var_level))
          ggsave(filename = paste0(output_dir_pubs, pathway, "_", contrast, "_model-", i, "_subset-var", subset_var, "_level-", subset_var_level, ".png"), height = 9, width = 10, units = "in")

        }
        
      } # End contrasts for loop.
    } # End model_numbers for loop.
  } # End subset_var_levels for loop.
} # End Subset variables for loop.

# Arrange plots into grid
# and add to PowerPoint.
# Graphing parameters.
plot_width <- 12
plot_height <- 12 
units <- "in"
res <- 300
plot_list_pathway_analysis_grid <- list()
for(sv in names(plot_list_pathway_analysis)) {
  for(svl in names(plot_list_pathway_analysis[[sv]])) {
    for(model_num in names(plot_list_pathway_analysis[[sv]][[svl]])) {
      # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
      # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
      # https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0

      p_list <- plot_list_pathway_analysis[[sv]][[svl]][[model_num]]

      n <- length(p_list)
      nCol <- ifelse(n <= 3, 2, floor(sqrt(n)))

      # Extract the legend from the plots.
      legend <- cowplot::get_plot_component(p_list[[1]], 'guide-box-right', return_all = TRUE)

      # Strip legends from p_list.
      for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}

      # Get the model information (test/contrast variable, random slope status, etc.)
      i <- model_number %>% regexPipes::gsub("model_", "")
      subset_var <- sv
      subset_var_level <- svl
      test_var <- results2 %>% dplyr::filter(`Model number`==i) %>% .$`Contrast variable` %>% .[1]

      if(subset_var=="DummySubsetVar") {
        subset_by <- ""
        filename_subset_by <- ""
      } else {
        subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
        filename_subset_by <- paste0("subset-variable-", subset_var, "_level-", subset_var_level)
      }

      # Arrange plots in p_list onto a grid.
      plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
      plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Pathway", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                                                         bottom = grid::textGrob("Normalized Enrichment Score", gp = grid::gpar(cex = 1.3)),
                                                         top = grid::textGrob(paste0("Pathway analysis ", subset_by, " \nTest (contrast) variable: ", test_var)))

      # Add back in the legend we extracted earlier.
      plot_grid2 <- grid.arrange(plot_grid, legend, ncol = 2, widths=c(10, 1))

      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, "FGSEA_bar-graphs_", filename_subset_by, "contrast-variable-", test_var, ".eps") #paste0(output_dir_pubs, "")
      png_path <- paste0(output_dir_pubs, "FGSEA_bar-graphs_", filename_subset_by, "contrast-variable-", test_var, ".png") #paste0(output_dir_pubs, "")
      plot <- plot_grid2 %>% ggpubr::as_ggplot()
      saveEPS(plot, eps_path, width = plot_width, height = plot_height)
      savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
     
      # ... add to the PowerPoint.
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Pathway analysis"),
                         location = ph_location_label(ph_label = "Title 1")) %>%
        officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                         location = ph_location_label(ph_label = "Content Placeholder 2"),
                         use_loc_size = FALSE)

      # Save to list.
      plot_list_pathway_analysis_grid[[sv]][[svl]][[model_num]] <- plot_grid2 %>% ggpubr::as_ggplot()

    }
  }
}

## ---------------------------
# Export to disk.

# Export PowerPoint file.
print(pptx, cl_args[5])
# Export NanoStringGeoMxSet as RDS file.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_pathway-analysis.rds"))
# Export FGSEA results as table.
pathway_df %>% write.csv(paste0(output_dir_tabular, "pathway-analysis_results.csv"))
# Export the raw plots as RDS file.
plot_list_pathway_analysis %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_raw-plots-list.rds"))
# Export the grid-arranged plots as RDS file.
plot_list_pathway_analysis_grid %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_grid-arranged-plots-list.rds"))