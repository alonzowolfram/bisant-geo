# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#     IMPORTANT NOTE ABOUT THIS MODULE:
#     For some weird reason, manipulation of the pData data frame in this module
#     messes up the rownames (including the sample names of the NanoStringGeoMx object).
#     I haven't figured out how to fix it yet, so for now we are just working around the issue. 
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
target_data_object_list <- readRDS(cl_args[5])
# We'll only need the main module for this one.
target_data_object <- target_data_object_list[[main_module]]

# Sanity check: number of compartment variables
print(compartment_vars)

# Update defaults for umap to contain a stable random_state (seed).
custom_umap <- umap::umap.defaults
custom_umap$random_state <- random_seed

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Unsupervised clustering ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Because of some weird thing with, I'm guessing, nested data frames,
# R thinks that the row names aren't unique whenever we manipulate pData.
# Right now, we're getting an issue when we try to reset the column names after unsupervised analysis.
# So we'll unnest the nested data frames (LOQ in our case.)
pData(target_data_object) <- pData(target_data_object) %>% unnest(cols = c(LOQ), names_sep = "-")
orig_var_names <- colnames(pData(target_data_object))
orig_row_names <- rownames(pData(target_data_object))

# Initialize the list to hold the plots.
plot_list_unsupervised_clustering <- list()
# Initialize the list to hold the dimension-reduction results.
dim_red_list <- list()

# For each normalization method, run UMAP, t-SNE, and PCA.
# https://www.biostars.org/p/9540137/
for(norm_method in names(target_data_object@assayData)) {
  if(norm_method %in% c("exprs", "bg_sub") | base::grepl("log_", norm_method)) next # Skip over raw data, background-subtracted-only data, and log-transformed data.

  message(paste0("Performing dimensionality reduction on data normalized with method ", norm_method))

  plot_list_unsupervised_clustering[[norm_method]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["UMAP"]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["t-SNE"]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["PCA"]] <- list()

  # Add 1 to expression matrix so we can log2 transform.
  exprs_mat <- assayDataElement(target_data_object, elt = norm_method) + 1
  # Check if any of the probes have 0 variance ... this will mess up PCA scaling if so.
  # https://stackoverflow.com/a/40317343/23532435
  # which(apply(exprs_mat, 1, var)==0)
  probes_0_var <- which(apply(exprs_mat, 1, var)==0)
  if(length(probes_0_var) > 0) {
    print("The following probes with 0 variance were detected and will be removed for dimensionality reduction: ")
    print(rownames(exprs_mat)[probes_0_var])
    exprs_mat <- exprs_mat[-probes_0_var,]
  }
  
  # Run UMAP.
  if(perform_UMAP) {
    umap_out <-
      umap(t(log2(exprs_mat)),
           config = custom_umap)
    umap_col_names <- paste0("UMAP", 1:2, "_", norm_method)
    pData(target_data_object)[, umap_col_names] <- umap_out$layout[, c(1,2)]
    dim_red_list[[paste0("UMAP_", norm_method)]] <- pData(target_data_object)[, umap_col_names]
  }
  
  # Run t-SNE.
  if(perform_tSNE) {
    set.seed(random_seed)
    tsne_out <-
      Rtsne(t(log2(exprs_mat)),
            perplexity = ncol(target_data_object)*.15)
    tsne_col_names <- paste0("tSNE", 1:2, "_", norm_method)
    pData(target_data_object)[, tsne_col_names] <- tsne_out$Y[, c(1,2)]
    dim_red_list[[paste0("tSNE_", norm_method)]] <- pData(target_data_object)[, tsne_col_names]
  }
  
  # Run PCA.
  # https://www.statology.org/principal-components-analysis-in-r/
  if(perform_PCA) {
    pca_out <- prcomp(t(log2(exprs_mat)), scale = TRUE)
    pca_out$rotation <- -1*pca_out$rotation
    pca_col_names <- paste0("PCA", 1:2, "_", norm_method)
    pData(target_data_object)[, pca_col_names] <- pca_out$x[, c(1,2)]
    dim_red_list[[paste0("PCA_", norm_method)]] <- pData(target_data_object)[, pca_col_names]
  }
  
  # For column renaming and removing later.
  n_reduction_types = perform_UMAP + perform_tSNE + perform_PCA
  n_cols_remove = 2 * n_reduction_types
  
  # Create a separate graph for each compartment variable.
  for(compartment_var in compartment_vars) {
    print(paste0("Compartment variable: ", compartment_var))

    # Check if the compartment_var is a "composite" variable (made up of 2+ variables).
    if(base::grepl("\\+", compartment_var)) {
      # Yes - so separate it into its component variables.
      compartment_var_comps <- compartment_var %>% stringr::str_split("\\+") %>% unlist
      compartment_var_1 <- compartment_var_comps[1]
      compartment_var_2 <- compartment_var_comps[2]

      # Rename the first compartment variable as "CompartmentVar1"
      # and the second as "CompartmentVar2."
      pData(target_data_object) <- pData(target_data_object) %>%
        dplyr::rename(CompartmentVar1 = !!as.name(compartment_var_1), CompartmentVar2 = !!as.name(compartment_var_2))
      pData(target_data_object)$CompartmentVar1 <- pData(target_data_object)$CompartmentVar1 %>% as.character %>% as.factor
      pData(target_data_object)$CompartmentVar2 <- pData(target_data_object)$CompartmentVar2 %>% as.character %>% as.factor

      # # Create a new column consisting of the component variables.
      # for(i in 1:length(compartment_var_comps)) {
      #   comp <- compartment_var_comps[i]
      #   if(i==1) {
      #     merged_col <- pData(target_data_object)[[comp]]
      #   } else {
      #     merged_col <- paste(merged_col, pData(target_data_object)[[comp]], sep = " | ")
      #   }
      # }
      # pData(target_data_object)[["CompartmentVar1"]] <- merged_col %>% as.character %>% as.factor
      # # But remember that we'll have to rename the columns to the original,
      # # so we'll have to remove this new column at the end before we rename.

    } else {
      # No - so rename the (single) compartment variable as "CompartmentVar1".
      pData(target_data_object) <- pData(target_data_object) %>%
        dplyr::rename(CompartmentVar1 = !!as.name(compartment_var))
      pData(target_data_object)$CompartmentVar1 <- pData(target_data_object)$CompartmentVar1 %>% as.character %>% as.factor
    }

    # UMAP
    if(perform_UMAP) {
      if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
        # CompartmentVar2 exists.
        plot_list_unsupervised_clustering[[norm_method]][["UMAP"]][[compartment_var]] <-
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(umap_col_names[1]), y = !!as.name(umap_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3, aes(shape = CompartmentVar2)) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var_1),
               shape = paste0(compartment_var_2)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      } else {
        # CompartmentVar2 doesn't exist.
        plot_list_unsupervised_clustering[[norm_method]][["UMAP"]][[compartment_var]] <-
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(umap_col_names[1]), y = !!as.name(umap_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      }
    }
    
    # t-SNE
    if(perform_tSNE) {
      if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
        plot_list_unsupervised_clustering[[norm_method]][["t-SNE"]][[compartment_var]] <-
          # CompartmentVar2 exists.
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(tsne_col_names[1]), y = !!as.name(tsne_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3, aes(shape = CompartmentVar2)) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var_1),
               shape = paste0(compartment_var_2)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      } else {
        plot_list_unsupervised_clustering[[norm_method]][["t-SNE"]][[compartment_var]] <-
          # CompartmentVar2 doesn't exist.
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(tsne_col_names[1]), y = !!as.name(tsne_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      }
    }
    
    # PCA
    if(perform_PCA) {
      if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
        plot_list_unsupervised_clustering[[norm_method]][["PCA"]][[compartment_var]] <-
          # CompartmentVar2 exists.
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(pca_col_names[1]), y = !!as.name(pca_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3, aes(shape = CompartmentVar2)) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var_1),
               shape = paste0(compartment_var_2)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      } else {
        plot_list_unsupervised_clustering[[norm_method]][["PCA"]][[compartment_var]] <-
          # CompartmentVar2 doesn't exist.
          ggplot(pData(target_data_object),
                 aes(x = !!as.name(pca_col_names[1]), y = !!as.name(pca_col_names[2]), color = CompartmentVar1)) +
          geom_point(size = 3) +
          labs(title = paste0("Compartment: ", compartment_var),
               color = paste0(compartment_var)) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank())
      }
    }
    
    # # If it was a composite variable, remove CompartmentVar before resetting the variable names.
    # if(ncol(pData(target_data_object)) > (length(orig_var_names) + 4)) { # + 4 to account for the t-SNE and UMAP cols.
    #   pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-CompartmentVar1)
    # }
    # Reset the variable names.
    colnames(pData(target_data_object))[1:(ncol(pData(target_data_object))-n_cols_remove)] <- c(orig_var_names)

    message(paste0("Success: compartment variable: ", compartment_var))
  } # End compartment variable loop.
  
  # Remove the last n columns (dimension reduction results)
  # where n = 2 * (perform_UMAP + perform_tSNE + perform_PCA)
  pData(target_data_object) <- pData(target_data_object) %>% .[,1:(ncol(.)-n_cols_remove)]
  
  message(paste0("Success: normalization method: ", norm_method))
} # End normalization method loop.

# Create the dimension reduction table and add to pData.
dim_red_df <- do.call(cbind, dim_red_list)
pData(target_data_object) <- cbind(pData(target_data_object), dim_red_df)

# # Export individual plots to EPS/PNG.
# plot_width = 8
# plot_height = 6
# units = "in"
# res = 300
# for(norm_method in names(plot_list_unsupervised_clustering)) {
#   for(dim_red_method in names(plot_list_unsupervised_clustering[[norm_method]])) {
#     for(compartment_var in names(plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]])) {
#       plot <- plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]][[compartment_var]]
# 
#       # # Set the scaling factors.
#       # # It's the same for both width and height since the grids are squares.
#       # scaling_factor <- ifelse(nCol > 1, (nCol / 2)^2, 1)  # Number of rows in current grid / 2 (base number)
# 
#       # Save to EPS and PNG and then ...
#       eps_path <- paste0(output_dir_pubs, 
#                          paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, "-", compartment_var, ".eps") %>% regexPipes::gsub("\\/", "_"))
#       png_path <- paste0(output_dir_pubs, 
#                          paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, "-", compartment_var, ".png") %>% regexPipes::gsub("\\/", "_"))
#       saveEPS(plot, eps_path, width = plot_width, height = plot_height)
#       savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
#     }
#   }
# }

# Arrange plots into grid.
plot_list_unsupervised_clustering_grid <- list()
for(norm_method in names(plot_list_unsupervised_clustering)) {
  message(paste0("Arranging plots for normalization method ", norm_method))
  
  plot_list_unsupervised_clustering_grid[[norm_method]] <- list()
  
  for(dim_red_method in names(plot_list_unsupervised_clustering[[norm_method]])) {
    if(length(plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]]) < 1) next
    
    message(paste0("Arranging plots for dimension reduction method ", dim_red_method, " for normalization method ", norm_method))
    
    p_list <- plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]]

    n <- length(p_list)
    nCol <- ifelse(n %in% 2:3, 2, floor(sqrt(n))) # If n = 1, floor(sqrt(n)) goes to 1.

    # Set the scaling factors for label and legend size.
    sqrt_n_col <- sqrt(nCol)
    scaling_factor <- ifelse(nCol > 1, (sqrt_n_col * nCol / 2), 1) # Number of rows in current grid / 2 (base number)

    # Arrange plots in p_list onto a grid.
    plot_grid <- do.call("grid.arrange", c(p_list, ncol=nCol))
    plot_grid <- plot_grid %>% ggpubr::annotate_figure(left = grid::textGrob("Significance, -log10(P)", hjust = 0, rot = 90, vjust = 1, gp = grid::gpar(cex = scaling_factor)),
                                                       bottom = grid::textGrob("", gp = grid::gpar(cex = scaling_factor)),
                                                       top = grid::textGrob(paste0(dim_red_method, " | Normalization: ", normalization_names[names(normalization_names)==norm_method]),
                                                                            gp = grid::gpar(cex = scaling_factor)
                                                                            )
                                                       )

    # Save to list.
    plot_list_unsupervised_clustering_grid[[norm_method]][[dim_red_method]] <- plot_grid
  }
}

# Add the module back to the list.
target_data_object_list[[main_module]] <- target_data_object

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## 16S scores ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Choose the dimension reductions to work with. 
dim_red_types <- c()
if(perform_UMAP) dim_red_types <- c(dim_red_types, "UMAP")
if(perform_tSNE) dim_red_types <- c(dim_red_types, "tSNE")
if(perform_PCA) dim_red_types <- c(dim_red_types, "PCA")

# Create a list to hold 16S score plots.
plot_list_16s_score <- list()

message("Creating graph of 16S scores.")
for(dim_red_type in dim_red_types) {
  
  # Map 16S scores onto dimension-reduction plot.
  for(norm_method in normalization_methods) {
    if(norm_method %in% c("exprs", "bg_sub") | base::grepl("log_", norm_method)) next # Skip over raw data, background-subtracted-only data, and log-transformed data.
    
    for(compartment_var in compartment_vars) {
      # Create the data frame.
      dim_red_col_names <- paste0(dim_red_type, 1:2, "_", norm_method)
      dat <- cbind(dim_red_list[[paste0(dim_red_type, "_", norm_method)]], 
        pData(target_data_object)[, c("Score16S", compartment_var)])
      
      # Create the graph.
      plot <- dat %>% 
        ggplot(aes(x = !!as.name(dim_red_col_names[1]), y = !!as.name(dim_red_col_names[2]))) + 
        geom_point(aes(size = Score16S, color = !!as.name(compartment_var))) +
        labs(title = paste0(dim_red_type, " | ", norm_method, " | ", "Compartment: ", compartment_var),
             color = paste0(compartment_var)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
      print(plot)
      
      plot_list_16s_score[[paste0(dim_red_type, "_", norm_method, "_", compartment_var)]] <- plot
      
    } # End compartment variable loop.
  } # End normalization methods loop.
} # End dimension reductions loop.

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Heatmap of high-CV genes ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# For each of the normalization methods, generate a heatmap.
plot_list_heatmap <- list()
cv_res <- c()
for(norm_method in names(target_data_object@assayData)) {
  if(norm_method %in% c("exprs", "bg_sub") | base::grepl("log_", norm_method)) next # Skip over raw data, background-subtracted-only data, and log-transformed data.

  print(paste0("Creating heatmap from data normalized with method ", norm_method))

  # Calculate coefficient of variance for each gene.
  CV_dat <- assayDataApply(target_data_object,
                           elt = paste0("log_", norm_method), MARGIN = 1, calc_CV)
  # Add to the data frame.
  cv_res <- cbind(cv_res, CV_dat)
  # Show the highest CD genes and their CV values.
  sort(CV_dat, decreasing = TRUE)[1:5]

  # Identify genes in the top 3rd of the CV values.
  GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8, na.rm = TRUE)] %>% .[!is.na(.)] # Remove NAs, probably the probes with 0 variance, like the negative controls.
  exprs_mat <- assayDataElement(target_data_object[GOI, ], elt = paste0("log_", norm_method))
  annot <- data.frame(pData(target_data_object)[, heatmap_ann_vars])
  rownames(annot) <- colnames(exprs_mat)
  plot_list_heatmap[[norm_method]] <- pheatmap(exprs_mat,
           scale = "row",
           show_rownames = FALSE, show_colnames = FALSE,
           border_color = NA,
           clustering_method = "average",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           breaks = seq(-3, 3, 0.05),
           color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
           annotation_col = annot, # https://www.researchgate.net/post/R-error-gpar-element-fill-must-not-be-length-0
           # Apparently pData loses the rownames, hence using the colnames of exprs_mat as the rownames of the annotation data frame.
           main = paste0("Normalization method: ", normalization_names[names(normalization_names)==norm_method])
           )

  message(paste0("Successfully created heatmap from data normalized with method ", norm_method))
}
colnames(cv_res) <- names(target_data_object@assayData) %>% .[!(. %in% c("exprs", "bg_sub") | base::grepl("log_", .))]

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Save the NanoStringGeoMxSet to RDS.
saveRDS(target_data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_unsupervised-analysis.rds"))
# Save the graphs to RDS.
saveRDS(plot_list_unsupervised_clustering, paste0(output_dir_rdata, "plot-list_unsupervised-clustering.rds"))
saveRDS(plot_list_unsupervised_clustering_grid, paste0(output_dir_rdata, "plot-list_unsupervised-clustering-grids.rds"))
saveRDS(plot_list_16s_score, paste0(output_dir_rdata, "plot-list_16s-score.rds"))
saveRDS(plot_list_heatmap, paste0(output_dir_rdata, "plot-list_cv_heatmaps.rds"))
# Save the CV results to CSV.
write.table(cv_res, paste0(output_dir_tabular, "CV_results_by-normalization-method.csv"), sep = ",")

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)