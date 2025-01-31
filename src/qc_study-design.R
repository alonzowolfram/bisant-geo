## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file.
source("src/setup.R")

# Read in the NanoStringGeoMxSet object. 
data_object_list <- readRDS(cl_args[5])

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Study design QC ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Access the PKC files, to ensure that expected PKCs have been loaded for this study.
modules <- names(data_object_list)
pkcs <- paste0(modules, ".pkc")
pkc_summary <- data.frame(PKCs = pkcs, modules = modules)

# # Visually summarize the experimental design for the dataset to look at the different types of samples and ROI/AOI segments; present in Sankey diagram.
# # Select the annotations we want to show, use `` to surround column names with
# # spaces or special symbols.
# count_mat <- count(pData(data_object), `slide name`, segment)
# # Gather the data and plot in order: class, slide name, region, segment.
# test_gr <- gather_set_data(count_mat, 1:4)
# test_gr$x <- factor(test_gr$x,
#                     levels = c("class", "slide name", "region", "segment"))
# # Plot Sankey.
# ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
#     geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
#     geom_parallel_sets_axes(axis.width = 0.2) +
#     geom_parallel_sets_labels(color = "white", size = 5) +
#     theme_classic(base_size = 17) +
#     theme(legend.position = "bottom",
#           axis.ticks.y = element_blank(),
#           axis.line = element_blank(),
#           axis.text.y = element_blank()) +
#     scale_y_continuous(expand = expansion(0)) +
#     scale_x_discrete(expand = expansion(0)) +
#     labs(x = "", y = "") +
#     annotate(geom = "segment", x = 4.25, xend = 4.25,
#              y = 20, yend = 120, lwd = 2) +
#     annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
#              hjust = 0.5, label = "100 segments")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save the NanoStringGeoMxSet object.
saveRDS(data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-study-design.rds"))
# Export the main module as a single object (not list) so we can use the Shiny app on it.
data_object <- data_object_list[[main_module]]
saveRDS(data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_raw_main-module.rds"))
# Save the PKC summary table.
saveRDS(pkc_summary, paste0(output_dir_rdata, "pkc_summary_table.rds"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)