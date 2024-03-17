## Source the setup.R file.
source("src/setup.R")

## ---------------------------
# Setup

# Read in the NanoStringGeoMxSet object. 
target_data_object <- readRDS(cl_args[4])
# Read in the PowerPoint.
pptx <- read_pptx(cl_args[5])

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
# Rename the compartment variable as "CompartmentVar".
pData(target_data_object) <- pData(target_data_object) %>% 
    dplyr::rename(CompartmentVar = !!as.name(compartment_var))

# Run UMAP.
# https://www.biostars.org/p/9540137/
umap_out <-
  umap(t(log2(assayDataElement(target_data_object, elt = normalization_method))),
       config = custom_umap)
pData(target_data_object)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
plot_list_unsupervised_clustering[["UMAP"]] <- ggplot(pData(target_data_object),
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
plot_list_unsupervised_clustering[["tSNE"]] <- ggplot(pData(target_data_object),
                                                      aes(x = tSNE1, y = tSNE2, color = CompartmentVar)) +
  geom_point(size = 3) +
  labs(title = paste0("t-SNE | Compartment = ", compartment_var)) +
  theme_bw()

# Reset the variable names.
colnames(pData(target_data_object)) <- c(orig_var_names, paste0("UMAP", 1:2), paste0("tSNE", 1:2))

# Save to RDS.
saveRDS(target_data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_unsupervised-analysis.rds"))

# Add a section header.
pptx <- pptx %>% 
  officer::add_slide(layout = "Section Header", master = "Office Theme") %>%
  officer::ph_with(value = paste0("Unsupervised clustering"), 
                   location = ph_location_label(ph_label = "Title 1"))
# Add the graphs.
for(item in names(plot_list_unsupervised_clustering)) {
  if(item != "neg_geo_means") {
    # Add to the PowerPoint. 
    pptx <- pptx %>%
      officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
      officer::ph_with(value = paste0("Unsupervised clustering"),
                       location = ph_location_label(ph_label = "Title 1")) %>% 
      officer::ph_with(value = plot_list_unsupervised_clustering[[item]],
                       location = ph_location_label(ph_label = "Content Placeholder 2"))
  } else {
    for(ann in names(plot_list_unsupervised_clustering[["neg_geo_means"]])) {
      # Add to the PowerPoint. 
      pptx <- pptx %>%
        officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
        officer::ph_with(value = paste0("Unsupervised clustering"),
                         location = ph_location_label(ph_label = "Title 1")) %>% 
        officer::ph_with(value = plot_list_unsupervised_clustering[[item]][[ann]],
                         location = ph_location_label(ph_label = "Content Placeholder 2"))
    }
  }
}
print(pptx, cl_args[5])