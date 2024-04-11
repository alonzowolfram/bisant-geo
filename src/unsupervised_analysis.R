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

# Initialize the list to hold the plots.
plot_list_unsupervised_clustering <- list()

# Update defaults for umap to contain a stable random_state (seed).
custom_umap <- umap::umap.defaults
custom_umap$random_state <- random_seed

# Because of some weird thing with, I'm guessing, nested data frames,
# R thinks that the row names aren't unique whenever we manipulate pData.
# Right now, we're getting an issue when we try to reset the column names after unsupervised analysis.
# So we'll unnest the nested data frames (LOQ in our case.)
pData(target_data_object) <- pData(target_data_object) %>% unnest(cols = c(LOQ), names_sep = "-")
orig_var_names <- colnames(pData(target_data_object))
orig_row_names <- rownames(pData(target_data_object))

# Run UMAP.
# https://www.biostars.org/p/9540137/
plot_list_unsupervised_clustering[["UMAP"]] <- list()
plot_list_unsupervised_clustering[["t-SNE"]] <- list()
for(compartment_var in compartment_vars) {
  # Check if the compartment_var is a "composite" variable (made up of 2+ variables).
  if(base::grepl("\\+", compartment_var)) {
    # Yes - so separate it into its component variables.
    compartment_var_comps <- compartment_var %>% stringr::str_split("\\+") %>% unlist
    # Create a new column consisting of the component variables.
    for(i in 1:length(compartment_var_comps)) {
      comp <- compartment_var_comps[i]
      if(i==1) {
        merged_col <- pData(target_data_object)[[comp]]
      } else {
        merged_col <- paste(merged_col, pData(target_data_object)[[comp]], sep = " | ")
      }
    }
    pData(target_data_object)[["CompartmentVar"]] <- merged_col
    # But remember that we'll have to rename the columns to the original, 
    # so we'll have to remove this new column at the end before we rename.

  } else {
    # No - so rename the (single) compartment variable as "CompartmentVar".
    pData(target_data_object) <- pData(target_data_object) %>% 
      dplyr::rename(CompartmentVar = !!as.name(compartment_var))
  }

  umap_out <-
  umap(t(log2(assayDataElement(target_data_object, elt = normalization_method))),
       config = custom_umap)
  pData(target_data_object)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
  plot_list_unsupervised_clustering[["UMAP"]][[compartment_var]] <- ggplot(pData(target_data_object),
                                                        aes(x = UMAP1, y = UMAP2, color = CompartmentVar)) +
                                                        geom_point(size = 3) +
                                                        labs(title = paste0("UMAP | Compartment = ", compartment_var)) +
                                                        theme_bw()

  # Run tSNE.
  set.seed(random_seed)
  tsne_out <-
    Rtsne(t(log2(assayDataElement(target_data_object, elt = normalization_method))),
          perplexity = ncol(target_data_object)*.15)
  pData(target_data_object)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
  plot_list_unsupervised_clustering[["tSNE"]][[compartment_var]] <- ggplot(pData(target_data_object),
                                                        aes(x = tSNE1, y = tSNE2, color = CompartmentVar)) +
    geom_point(size = 3) +
    labs(title = paste0("t-SNE | Compartment = ", compartment_var)) +
    theme_bw()

  # If it was a composite variable, remove CompartmentVar before resetting the variable names.
  if(ncol(pData(target_data_object)) > (length(orig_var_names) + 4)) { # + 4 to account for the t-SNE and UMAP cols.
    pData(target_data_object) <- pData(target_data_object) %>% dplyr::select(-CompartmentVar)
  }
  # Reset the variable names.
  colnames(pData(target_data_object)) <- c(orig_var_names, paste0("UMAP", 1:2), paste0("tSNE", 1:2))
}

# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Unsupervised clustering"), 
                   location = ph_location_label(ph_label = "Title 1"))

# Graphing parameters.
plot_width <- 8
plot_height <- 6
units <- "in"
res <- 300
# Add the graphs.
for(item in names(plot_list_unsupervised_clustering)) {
  for(item2 in names(plot_list_unsupervised_clustering[[item]])) {
    plot <- plot_list_unsupervised_clustering[[item]][[item2]]
    
    # Save to EPS and PNG and then ...
    eps_path <- paste0(output_dir_pubs, "unsupervised-analysis_", item, "-", item2, ".eps")
    png_path <- paste0(output_dir_pubs, "unsupervised-analysis_", item, "-", item2, ".png")
    saveEPS(plot, eps_path, width = plot_width, height = plot_height)
    savePNG(plot, png_path, width = plot_width, height = plot_height, units = units, res = res)
    
    # Add to the PowerPoint. 
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Unsupervised clustering"),
                       location = ph_location_label(ph_label = "Title 1")) %>% 
      officer::ph_with(value = external_img(png_path, width = plot_width, height = plot_height, unit = units),
                       location = ph_location_label(ph_label = "Content Placeholder 2"),
                       use_loc_size = FALSE)
    
  }
}

## ---------------------------
# Export to disk.

# Save the NanoStringGeoMxSet to RDS.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_unsupervised-analysis.rds"))
# Output everything to the PowerPoint. 
print(pptx, cl_args[5])