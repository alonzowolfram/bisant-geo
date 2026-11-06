## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
data_object_list <- readRDS(cl_args[5])

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Study design QC ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# PKC summary
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Access the PKC files, to ensure that expected PKCs have been loaded for this study
modules <- names(data_object_list)
pkcs <- paste0(modules, ".pkc")
# Create summary table
pkc_summary <- data.frame(PKCs = pkcs, modules = modules)

# Set `main_module` if not set already
if(flagVariable(main_module)) main_module <- modules[1]
# Set the data object
data_object <- data_object_list[[main_module]]

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Subject-level summary
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Create summary table at the individual subject level by specified (categorical) characteristics
subject_summary_tables <- list()

# If `subject_categories` is set, split it into its components by semicolons
if(!flagVariable(subject_categories)) {
  # Split `subject_categories` by semicolons
  # This will give groups of subject categories, each group corresponding to one entry in `subject_ids`
  subject_categories <- subject_categories %>% str_split(";") %>% unlist
}
# If `subject_ids` and `subject_categories` are set, create table for each individual subject category, for each subject identifier
if(!flagVariable(subject_ids) & !flagVariable(subject_categories)) {
  # Extract the pData data frame
  pdata <- pData(data_object)
  
  # Create a table of subject identifiers - subject categories
  # Recycle entries if one is longer than the other
  tables_df <- cbind(subject_ids, subject_categories)
  
  for(i in 1:nrow(tables_df)) {
    subject_id <- tables_df[i,1]
    subject_category <- tables_df[i,2] %>% unlist
    
    # For each entry in `subject_categories`, add to a _list_
    # splitting the entry into separate entries by commas (",") first and then slashes ("/")
    subject_categories_list <- subject_category %>% str_split(",") %>% .[[1]] %>%
                                  lapply(FUN = function(x) {
                                        x %>% str_split("/") %>% .[[1]]
                                      })
    
    for(entry in subject_categories_list) {
      # Create the identifier for this entry
      entry_cleaned <- paste(entry, collapse = "_")
      name <- glue::glue("{subject_id} by {entry_cleaned}")
      
      # Use the `pdata` data frame to create confusion matrix for each individual subject category
      subject_summary_tables[[name]] <- pdata %>% 
        dplyr::select(!!as.name(subject_id), entry) %>% 
        dplyr::distinct() %>% 
        dplyr::select(entry) %>% 
        table
    }
    
  }
  
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Save the NanoStringGeoMxSet object
saveRDS(data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_qc-study-design.rds"))
# Export the main module as a single object (not list) so we can use the Shiny app on it
data_object <- data_object_list[[main_module]]
saveRDS(data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_raw_main-module.rds"))
# Save the PKC summary table
saveRDS(pkc_summary, paste0(output_dir_rdata, "pkc_summary_table.rds"))
# Save the subject characteristics summary tables
saveRDS(subject_summary_tables, paste0(output_dir_rdata, "subject_summary_tables.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_qc_study-design.RData"))

# Update latest module completed
updateLatestModule(output_dir_rdata, current_module)