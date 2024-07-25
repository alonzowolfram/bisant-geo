## ---------------------------
#     IMPORTANT NOTE ABOUT THIS MODULE:
#     For some weird reason, manipulation of the pData data frame in this module
#     messes up the rownames (including the sample names of the NanoStringGeoMx object).
#     I haven't figured out how to fix it yet, so for now we are just working around the issue. 
## ---------------------------

## Source the setup.R file.
source("src/setup.R")

## ---------------------------
# Setup

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

# sanity check: number of compartment variables
print(compartment_vars)

# Update defaults for umap to contain a stable random_state (seed).
custom_umap <- umap::umap.defaults
custom_umap$random_state <- random_seed

# Add a section header.
pptx <- pptx %>%
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Unsupervised clustering"),
                   location = ph_location_label(ph_label = "Title 1"))

## ---------------------------
# Unsupervised clustering

# Because of some weird thing with, I'm guessing, nested data frames,
# R thinks that the row names aren't unique whenever we manipulate pData.
# Right now, we're getting an issue when we try to reset the column names after unsupervised analysis.
# So we'll unnest the nested data frames (LOQ in our case.)
pData(target_data_object) <- pData(target_data_object) %>% unnest(cols = c(LOQ), names_sep = "-")
orig_var_names <- colnames(pData(target_data_object))
orig_row_names <- rownames(pData(target_data_object))

# Initialize the list to hold the plots.
plot_list_unsupervised_clustering <- list()

# For each normalization method, run UMAP, t-SNE, and PCA.
# https://www.biostars.org/p/9540137/
for(norm_method in names(target_data_object@assayData)) {
  if(norm_method %in% c("exprs", "bg_sub", "log_norm")) next # Skip over raw data and background-subtracted-only data.

  print(paste0("Performing dimensionality reduction on data normalized with method ", norm_method))

  plot_list_unsupervised_clustering[[norm_method]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["UMAP"]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["t-SNE"]] <- list()
  plot_list_unsupervised_clustering[[norm_method]][["PCA"]] <- list()

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
    umap_out <-
      umap(t(log2(exprs_mat)),
           config = custom_umap)
    pData(target_data_object)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
    if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
      plot_list_unsupervised_clustering[[norm_method]][["UMAP"]][[compartment_var]] <-
        # CompartmentVar2 exists.
        ggplot(pData(target_data_object),
               aes(x = UMAP1, y = UMAP2, color = CompartmentVar1)) +
        geom_point(size = 3, aes(shape = CompartmentVar2)) +
        labs(title = paste0("Compartment: ", compartment_var),
             color = paste0(compartment_var_1),
             shape = paste0(compartment_var_2)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    } else {
      plot_list_unsupervised_clustering[[norm_method]][["UMAP"]][[compartment_var]] <-
        # CompartmentVar2 doesn't exist.
        ggplot(pData(target_data_object),
               aes(x = UMAP1, y = UMAP2, color = CompartmentVar1)) +
        geom_point(size = 3) +
        labs(title = paste0("Compartment: ", compartment_var),
             color = paste0(compartment_var)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    }

    # Run t-SNE.
    set.seed(random_seed)
    tsne_out <-
      Rtsne(t(log2(exprs_mat)),
            perplexity = ncol(target_data_object)*.15)
    pData(target_data_object)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
    if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
      plot_list_unsupervised_clustering[[norm_method]][["t-SNE"]][[compartment_var]] <-
        # CompartmentVar2 exists.
        ggplot(pData(target_data_object),
               aes(x = tSNE1, y = tSNE2, color = CompartmentVar1)) +
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
               aes(x = tSNE1, y = tSNE2, color = CompartmentVar1)) +
        geom_point(size = 3) +
        labs(title = paste0("Compartment: ", compartment_var),
             color = paste0(compartment_var)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    }

    # Run PCA.
    # https://www.statology.org/principal-components-analysis-in-r/
    pca_out <- prcomp(t(log2(exprs_mat)), scale = TRUE)
    pca_out$rotation <- -1*pca_out$rotation
    pData(target_data_object)[, c("PCA1", "PCA2")] <- pca_out$x[, c(1,2)]
    if("CompartmentVar2" %in% colnames(pData(target_data_object))) {
      plot_list_unsupervised_clustering[[norm_method]][["PCA"]][[compartment_var]] <-
        # CompartmentVar2 exists.
        ggplot(pData(target_data_object),
               aes(x = PCA1, y = PCA2, color = CompartmentVar1)) +
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
               aes(x = PCA1, y = PCA2, color = CompartmentVar1)) +
        geom_point(size = 3) +
        labs(title = paste0("Compartment: ", compartment_var),
             color = paste0(compartment_var)) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    }

    # # If it was a composite variable, remove CompartmentVar before resetting the variable names.
    # if(ncol(pData(target_data_object)) > (length(orig_var_names) + 4)) { # + 4 to account for the t-SNE and UMAP cols.
    #   pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-CompartmentVar1)
    # }
    # Reset the variable names.
    colnames(pData(target_data_object)) <- c(orig_var_names, paste0("UMAP", 1:2), paste0("t-SNE", 1:2), paste0("PCA", 1:2))
    # Remove the last 6 columns (dimension reduction results).
    pData(target_data_object) <- pData(target_data_object) %>% .[,1:(ncol(.)-6)]

    print(paste0("Success: compartment variable: ", compartment_var))
  }
  print(paste0("Success: normalization method: ", norm_method))
}

# Export individual plots to EPS/PNG.
plot_width = 8
plot_height = 6
units = "in"
res = 300
for(norm_method in names(plot_list_unsupervised_clustering)) {
  for(dim_red_method in names(plot_list_unsupervised_clustering[[norm_method]])) {
    for(compartment_var in names(plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]])) {
      plot <- plot_list_unsupervised_clustering[[norm_method]][[dim_red_method]][[compartment_var]]

      # # Set the scaling factors.
      # # It's the same for both width and height since the grids are squares.
      # scaling_factor <- ifelse(nCol > 1, (nCol / 2)^2, 1)  # Number of rows in current grid / 2 (base number)

      # Save to EPS and PNG and then ...
      eps_path <- paste0(output_dir_pubs, 
                         paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, "-", compartment_var, ".eps") %>% regexPipes::gsub("\\/", "_"))
      png_path <- paste0(output_dir_pubs, 
                         paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, "-", compartment_var, ".png") %>% regexPipes::gsub("\\/", "_"))
      saveEPS(plot, eps_path, width = plot_width, height = plot_height)
      savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
    }
  }
}

# Arrange plots into grid.
plot_list_unsupervised_clustering_grid <- list()
for(norm_method in names(plot_list_unsupervised_clustering)) {
  for(dim_red_method in names(plot_list_unsupervised_clustering[[norm_method]])) {
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
    plot_list_unsupervised_clustering_grid[[norm_method]][[dim_red_method]][["plot"]] <- plot_grid
    plot_list_unsupervised_clustering_grid[[norm_method]][[dim_red_method]][["nCol"]] <- nCol
  }
}

# Graphing parameters.
plot_width <- 18
plot_height <- 20
units <- "in"
res <- 300
# Add the dimensionality reduction (PCA, UMAP, t-SNE) graphs.
for(norm_method in names(plot_list_unsupervised_clustering)) {
  for(dim_red_method in names(plot_list_unsupervised_clustering[[norm_method]])) {
    plot <- plot_list_unsupervised_clustering_grid[[norm_method]][[dim_red_method]][["plot"]]
    nCol <- plot_list_unsupervised_clustering_grid[[norm_method]][[dim_red_method]][["nCol"]]

    # Set the scaling factors.
    # It's the same for both width and height since the grids are squares.
    scaling_factor <- ifelse(nCol > 1, (nCol / 2)^2, 1)  # Number of rows in current grid / 2 (base number)

    # Save to EPS and PNG and then ...
    eps_path <- paste0(output_dir_pubs, 
                       paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, ".eps") %>% regexPipes::gsub("\\/", "_"))
    png_path <- paste0(output_dir_pubs, 
                       paste0("unsupervised-analysis_", norm_method, "-", dim_red_method, ".png") %>% regexPipes::gsub("\\/", "_"))
    saveEPS(plot, eps_path, width = plot_width, height = plot_height)
    savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)

    # Add to the PowerPoint.
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Unsupervised clustering | ", normalization_names[names(normalization_names)==norm_method], " | ", dim_red_method),
                       location = ph_location_label(ph_label = "Title 1")) %>%
      officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                       location = ph_location_label(ph_label = "Content Placeholder 2"),
                       use_loc_size = FALSE)
  }
}

## ---------------------------
# Heatmap of high-CV genes

# For each of the normalization methods, generate a heatmap.
plot_list_heatmap <- list()
cv_res <- c()
for(norm_method in names(target_data_object@assayData)) {
  if(norm_method %in% c("exprs", "bg_sub", "log_norm")) next # Skip over raw data, background-subtracted-only data, and the log-norm data from the previous iteration.

  print(paste0("Creating heatmap from data normalized with method ", norm_method))

  # Create a log2 transform of the data for analysis.
  assayDataElement(object = target_data_object, elt = "log_norm", validate = FALSE) <- # Have to set validate to FALSE; otherwise it thinks the dimensions aren't the same. ...
    assayDataApply(target_data_object, 2, FUN = function(x) log2(x+1), elt = norm_method)

  # Calculate coefficient of variance for each gene.
  CV_dat <- assayDataApply(target_data_object,
                           elt = "log_norm", MARGIN = 1, calc_CV)
  # Add to the data frame.
  cv_res <- cbind(cv_res, CV_dat)
  # Show the highest CD genes and their CV values.
  sort(CV_dat, decreasing = TRUE)[1:5]

  # Identify genes in the top 3rd of the CV values.
  GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8, na.rm = TRUE)] %>% .[!is.na(.)] # Remove NAs, probably the probes with 0 variance, like the negative controls.
  exprs_mat <- assayDataElement(target_data_object[GOI, ], elt = "log_norm")
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

  print(paste0("Successfully created heatmap from data normalized with method ", norm_method))
}
colnames(cv_res) <- names(target_data_object@assayData) %>% .[!(. %in% c("exprs", "bg_sub", "log_norm"))]

## ---------------------------
# Output

# Graphing parameters.
plot_width <- 8
plot_height <- 8
units <- "in"
res <- 300

# Add the heatmaps.
for(norm_method in names(target_data_object@assayData)) {
  if(norm_method %in% c("exprs", "bg_sub", "log_norm")) next # Skip over raw data, background-subtracted-only data, and the log-norm data from the previous iteration.

  plot <- plot_list_heatmap[[norm_method]]

  # Save to EPS and PNG and then ...
  eps_path <- paste0(output_dir_pubs, "heatmap_high-CV-genes_", norm_method, "-normalized.eps")
  png_path <- paste0(output_dir_pubs, "heatmap_high-CV-genes_", norm_method, "-normalized.png")
  saveEPS(plot, eps_path, width = plot_width, height = plot_height)
  savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)

  # Add to the PowerPoint.
  pptx <- pptx %>%
    officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
    officer::ph_with(value = paste0("Heatmap of high-CV genes"),
                     location = ph_location_label(ph_label = "Title 1")) %>%
    officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                     location = ph_location_label(ph_label = "Content Placeholder 2"),
                     use_loc_size = FALSE)
}

## ---------------------------
# Export to disk.

# Output everything to the PowerPoint.
print(pptx, cl_args[5])
# Save the NanoStringGeoMxSet to RDS.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_unsupervised-analysis.rds"))
# Save the graphs to RDS.
saveRDS(plot_list_unsupervised_clustering, paste0(output_dir_rdata, "plot-list_unsupervised-clustering.rds"))
saveRDS(plot_list_unsupervised_clustering_grid, paste0(output_dir_rdata, "plot-list_unsupervised-clustering-grids.rds"))
saveRDS(plot_list_heatmap, paste0(output_dir_rdata, "plot-list_heatmaps.rds"))
# Save the CV results to CSV.
write.table(cv_res, paste0(output_dir_tabular, "CV_results_by-normalization-method.csv"), sep = ",")