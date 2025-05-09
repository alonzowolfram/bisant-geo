## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Setting up for pathway analysis.")

## Source the setup.R file.
source("src/setup.R")

# Read in the DE genes table.
message("Reading in DE genes table.")
results2 <- read.csv(cl_args[5], row.names = 1, check.names = FALSE)
print(paste0("results2 dimensions: ", dim(results2)))

# If pathway_table_file is provided, check that it's a valid file.
# If not, default to hallmark, BioCarta, and Reactome pathways.
# In future iterations, we'll probably include a default pathway_table_file instead of hard-coding the pathways here.
set.default.pathways <- function() {
  cats <- c("H") # , rep("C2", 2)
  subcats <- c(NA) # , "CP:BIOCARTA", "CP:REACTOME"
  return(list(cats = cats, subcats = subcats))
}

if(!flagVariable(pathway_table_file)) { 
  if(!file.exists(pathway_table_file)) {
    msigdb_list <- set.default.pathways()
  } else {
    # Read in the file and check that it has at least one entry. 
    message("Checking provided pathway file.")
    if(base::grepl("\\.csv$", pathway_table_file)) { pathway_table <- read.csv(pathway_table_file, header = F)}
    else if(base::grepl("\\.tsv$", pathway_table_file)) { pathway_table <- read.table(pathway_table_file, header = F, sep = "\t") }
    else if(base::grepl("\\.xls.*$", pathway_table_file)) { pathway_table <- read_excel(pathway_table_file, col_names = F, na = "NA")} 
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

# Set graphical parameters.
nes_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Set the normalization method.
normalization_method <- normalization_names[names(normalization_names)==normalization_methods[1]]

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Pathway analysis (FGSEA) ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## ................................................
##
### Data preparation ----
##
## ................................................
# Subset results2 table (data frame of differentially expressed genes).
message(paste0("Subsetting results to include only those from the normalization method ", normalization_method, "."))
results2_sub <- results2 %>% dplyr::filter(`Normalization method`==normalization_method)

# Load the gene sets (from MSigDB).
# https://rpubs.com/LiYumei/806213
message("Loading gene sets.")
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
message("Subsetting gene sets.")
if(sum(!is.null(individual_pathways)) > 0 & sum(individual_pathways != "") > 0) gene_sets <- gene_sets %>% dplyr::filter(gs_name %in% individual_pathways)

# Create an MSigDB pathway list.
message("Creating MSigDB gene set list.")
msigdbr_pathway_list <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)

## ................................................
##
### Calculation and graphing ----
##
## ................................................
# Loop through DE gene table and perform FGSEA.
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
plot_list_pathway_analysis <- list()
pathway_df <- data.frame()
models <- results2_sub$`Model` %>% unique
model_number <- 1
model_list <- list() # Hold all the model formulae so we can access them when arranging the plots into grids.

message("Performing FGSEA.")
for(subset_var in unique(results2_sub$`Subset variable`)) { # We're not naming it subset_vars because we already have a variable by that name. ... 
  
  subset_var_levels <- results2_sub %>% dplyr::filter(`Subset variable`==subset_var) %>% .$`Subset level` %>% unique
  for(subset_var_level in subset_var_levels) {
    
    for(model in models) {
      # Get the contrasts.
      contrasts <- results2_sub %>% dplyr::filter(`Subset variable`==subset_var & `Subset level`==subset_var_level & `Model`==model) %>% .$Contrast %>% unique
      # Get the contrast variable. 
      contrast_var <- results2_sub %>% dplyr::filter(`Model`==model) %>% .$`Contrast variable` %>% .[1]
      
      for(contrast in contrasts) {
        message(paste0("Subset variable ", subset_var, " | level ", subset_var_level, " | model ", model, " | contrast ", contrast))
        # Subset. 
        results2_sub_sub <- results2_sub %>% dplyr::filter(
          `Subset variable`==subset_var &
            `Subset level`==subset_var_level &
            `Model`==model &
            Contrast==contrast
        )
        
        # Extract the DEG data--these will be the ranks we use as input to FGSEA. 
        ranks <- results2_sub_sub$Estimate
        names(ranks) <- results2_sub_sub$Gene
        ranks <- ranks %>% .[order(.)]
        
        #### Calculation ----
        message("Running FGSEA.")
        fgseaRes <- fgsea(pathways = msigdbr_pathway_list, 
                          stats    = ranks,
                          minSize  = 1, # 15
                          maxSize  = Inf) # 500
        df_sub <- fgseaRes %>% dplyr::filter(is.finite(NES))
        # Make the pathway names more readable.
        message("Cleaning FGSEA results.")
        df_sub$PathwayCleaned <- df_sub$pathway %>% 
          stringr::str_split("_") %>% 
          lapply(FUN = function(x) c(paste0(x[1], ":"), x[2:length(x)]) %>% paste(collapse=" ")) %>%  
          unlist() %>% 
          stringr::str_split(": ") %>% 
          lapply(FUN = function(x) stringr::str_to_sentence(x) %>% paste(collapse=": ")) %>% 
          unlist()
        
        # Add pathway ranking score (-log10(padj) * |NES|)
        message("Adding pathway ranking scores.")
        
        df_sub$PathwayScore <- -log10(df_sub$padj) * abs(df_sub$NES)
        # Before calculating percentiles for pathway, cull pathway list down to its final form (i.e., what will actually be graphed.)
        # If individual_pathways is set, subset to include only the pathways of interest.
        if(!is.null(individual_pathways) & sum(individual_pathways=="") < length(individual_pathways)) df_sub <- df_sub %>% dplyr::filter(pathway %in% individual_pathways)
        
        # Add information about the model (model number, contrast variable, current contrast.)
        message("Adding information about model.")
        df_sub$`Model` <- model
        df_sub$`Contrast variable` <- contrast_var
        df_sub$Contrast <- contrast
          
        # Add column indicating percentile (how high up the list a given pathway is.)
        df_sub <- df_sub %>% 
          dplyr::arrange(-PathwayScore) %>% 
          dplyr::mutate(percentile = 100 - (1:nrow(.) / nrow(.) * 100))
        df_sub$percentile <- format(round(df_sub$percentile, 1), nsmall = 1)
        
        ### Graphing ----
        message("Graphing.")
        
        upper_limit <- ifelse(max(df_sub$NES) > 0, max(df_sub$NES) %>% ceiling(), max(df_sub$NES) %>% abs %>% ceiling())
        lower_limit <- ifelse(min(df_sub$NES) < 0, min(df_sub$NES) %>% floor(), min(df_sub$NES) %>% -. %>% floor())
        upper_limit_axis <- ifelse(max(df_sub$NES) < 0, 0, upper_limit)
        lower_limit_axis <- ifelse(min(df_sub$NES) > 0, 0, lower_limit)
        dodgewidth <- position_dodge(width=0.9)
        
        # Get the contrast elements for the plot label.
        contrast_element_1 <- contrast %>% strsplit(" - ") %>% unlist %>% .[1]
        contrast_element_2 <- contrast %>% strsplit(" - ") %>% unlist %>% .[2]
        
        # Graphing parameters.
        # Bar width needs to _increase_ with the # of rows in df_sub. 0.75 is a good width for 15 rows, and 0.25 is a good width for <= 3 rows, so ... 
        # sqrt function? Or sigmoidal?
        # We'll go with a piecewise function.
        # Code base stolen from https://stackoverflow.com/a/8788595/23532435
        # Second piece of the function will be linear (y = mx + b).
        # Slope for second piece of the function: rise / run = (0.75 - 0.25) / (5 - 1) = 0.042
        # y-intercept: b = y - mx
        slope <- (0.75 - 0.25) / (5 - 1)
        y_int <- 0.75 - slope * 5
        ratio <- nrow(df_sub) / 3
        bar_width <- (ratio <= 1) * 0.25 + 
                     (ratio > 1 & ratio <= 5) * (slope * ratio + y_int) + 
                     (ratio > 5) * ((sqrt(ratio) - sqrt(5)) + 0.75)
        
        # Get the data frame into the proper format for graphing using geom_rect.
        df_sub_graphing <- df_sub 
        
        df_sub_graphing$PathwayCleaned <- df_sub_graphing$PathwayCleaned %>% as.factor %>% factor(levels = df_sub %>% dplyr::arrange(NES) %>% .$PathwayCleaned) # Arrange factor levels of pathway by NES.
        df_sub_graphing$ymin <- ifelse(df_sub_graphing$NES >= 0, 0, df_sub_graphing$NES)
        df_sub_graphing$ymax <- ifelse(df_sub_graphing$NES >= 0, df_sub_graphing$NES, 0)
        df_sub_graphing$xmin <- (df_sub_graphing$PathwayCleaned %>% as.numeric) - (bar_width/2)
        df_sub_graphing$xmax <- (df_sub_graphing$PathwayCleaned %>% as.numeric) + (bar_width/2)
        
        # If the number of pathways > limits for graphing, perform cutoff by PathwayScore.
        if(!is.null(n_max_pathways) & n_max_pathways != "") if(nrow(df_sub_graphing) > n_max_pathways) df_sub_graphing <- df_sub_graphing %>% dplyr::top_n(n = n_max_pathways, wt = PathwayScore)
        
        # Truncate pathway names (doing it here, AFTER conversion to factor, because if truncating makes two pathway names come out the same, we can't have duplicated factors)
        # Reorder factor levels by NES (so the bars are sorted)
        df_sub_graphing$PathwayCleanedChar <- df_sub_graphing$PathwayCleaned %>%
          as.character() %>%
          truncate_strings()
        
        df_sub_graphing$PathwayCleanedChar <- factor(
          df_sub_graphing$PathwayCleanedChar,
          levels = df_sub_graphing %>% arrange(NES) %>% pull(PathwayCleanedChar)
        )
        
        # Graph.
        plot <- ggplot(df_sub_graphing, aes(x = PathwayCleanedChar, y = NES, fill = NES)) +
          geom_col(width = 0.8) +  # Adjust width as needed
          scale_fill_gradientn(colors = nes_palette, limits = c(lower_limit, upper_limit)) +
          coord_flip() +
          theme_bw() +
          theme(
            aspect.ratio = 1/1,
            legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks.y = element_blank()
          ) +
          labs(x = "", y = paste0(contrast_element_2, " <-> ", contrast_element_1))
        
        # plot <- ggplot(df_sub_graphing, aes(y = NES, x = reorder(PathwayCleanedChar, NES), fill = NES)) +
        #   # Bar graph manually, using geom_rect(angle).
        #   geom_rect(xmin = df_sub_graphing$xmin, xmax = df_sub_graphing$xmax, ymin = df_sub_graphing$ymin, ymax = df_sub_graphing$ymax) +
        #   # Set the color scale for the bars.
        #   scale_fill_gradientn(colors = nes_palette, limits = c(lower_limit, upper_limit)) + 
        #   # Flip horizontal.
        #   coord_flip() + 
        #   # B&W theme.
        #   theme_bw() + 
        #   # Remove grid lines and set aspect ratio.
        #   theme(aspect.ratio = 1/1, # Aspect ratio needs to _decrease_ (limit 0) with the # of rows in df_sub_graphing.
        #         legend.position = "right",
        #         panel.grid.minor = element_blank(),
        #         panel.grid.major = element_blank(),
        #         # axis.text.x=element_blank(), 
        #         # axis.ticks.x=element_blank(), 
        #         # axis.text.y=element_blank(), 
        #         axis.ticks.y=element_blank()) + 
        #   # Set the labels.
        #   labs(x = "", y = paste0(contrast_element_2, " <-> ", contrast_element_1))
        
        # Add to the list.
        plot_list_pathway_analysis[[subset_var]][[subset_var_level]][[paste0("model_", model_number)]][[contrast]] <- plot
        
        # Add df_sub to pathway_df
        pathway_df <- rbind(pathway_df, df_sub %>% dplyr::select(-leadingEdge))

        # Graph enrichment plots.
        message("Graphing enrichment plots.")
        for(pathway in names(msigdbr_pathway_list)) {
          plotEnrichment(msigdbr_pathway_list[[pathway]], ranks, gseaParam = 1, ticksSize = 0.2) +
            ggtitle(paste0(pathway, " | contrast: ", contrast, " | subset by ", subset_var, " | level: ", subset_var_level))
          ggsave(filename = paste0(output_dir_pubs, 
                                   paste0(pathway, "_", contrast, "_model-", model_number, "_subset-var", subset_var, "_level-", subset_var_level, ".png") %>% regexPipes::gsub("\\/", "_")), 
                 height = 9, width = 10, units = "in")

        }
        
      } # End contrasts for loop.
      
      model_list[[paste0("model_", model_number)]] <- as.character(model)
      model_number <- model_number + 1
      
    } # End models for loop.
  } # End subset_var_levels for loop.
} # End Subset variables for loop.

## ................................................
##
### Arranging plots into grids ----
##
## ................................................
# Graphing parameters.
plot_width <- 20
plot_height <- 15
units <- "in"
res <- 300

plot_list_pathway_analysis_grid <- list()
message("Arranging plots into grid.")
for(sv in names(plot_list_pathway_analysis)) {
  for(svl in names(plot_list_pathway_analysis[[sv]])) {
    for(model_num in names(plot_list_pathway_analysis[[sv]][[svl]])) {
      # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
      # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
      # https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0

      p_list <- plot_list_pathway_analysis[[sv]][[svl]][[model_num]]

      n <- length(p_list)
      if(n > 16) {
        # Stop right there.
        warning(paste0("The combination of ", sv, " - ", svl, " - ", model_num, " has ", n, " plots, too many for graphing. Please graph these manually. Skipping to the next list of graphs."))
        rm(p_list)
        gc()
        next
      }
      nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.
      
      # Set the scaling factors for label and legend size.
      sqrt_n_col <- sqrt(nCol)
      scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 3), 0.35) # Number of rows in current grid / 3 (base number)
      res_scaling_factor <- max(scaling_factor, 1)

      # Each plot will have a different scale, so we will not include a common legend.
      # for(item in names(p_list)) {p_list[[item]] <- p_list[[item]] + theme(legend.position = "none")}

      # Get the model information (test/contrast variable, etc.)
      model <- model_list[[model_num]]
      subset_var <- sv
      subset_var_level <- svl
      test_var <- results2_sub %>% dplyr::filter(`Model`==model) %>% .$`Contrast variable` %>% .[1]

      if(subset_var=="DummySubsetVar") {
        subset_by <- ""
        filename_subset_by <- ""
      } else {
        subset_by <- paste0("| Subset variable: ", subset_var, ", level: ", subset_var_level)
        filename_subset_by <- paste0("subset-variable-", subset_var, "_level-", subset_var_level, "_")
      }

      # Arrange plots in p_list onto a grid.
      plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
      plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Pathway", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                                                         bottom = grid::textGrob("Normalized Enrichment Score", gp = grid::gpar(cex = 1.3)),
                                                         top = grid::textGrob(paste0("Pathway analysis ", subset_by, 
                                                                                     "\nTest (contrast) variable: ", test_var,
                                                                                     "\nNormalization method: ", normalization_method)))

      # Save to list.
      plot_list_pathway_analysis_grid[[sv]][[svl]][[model_num]] <- plot_grid#2 %>% ggpubr::as_ggplot()

    }
  }
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message("Exporting to disk.")

# Export FGSEA results as table.
pathway_df %>% write.csv(paste0(output_dir_tabular, "pathway-analysis_results.csv"))
# Export the raw plots as RDS file.
plot_list_pathway_analysis %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_raw-plots-list.rds"))
# Export the grid-arranged plots as RDS file.
plot_list_pathway_analysis_grid %>% saveRDS(paste0(output_dir_rdata, "pathway-analysis_grid-arranged-plots-list.rds"))