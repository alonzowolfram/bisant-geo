## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file.
source("src/setup.R")

# Automatically list files in each directory for use.
dcc_files <- dir(dcc_dir, pattern = ".dcc$",
                 full.names = TRUE, recursive = TRUE)
if(!is.null(pkc_filenames) && sum(pkc_filenames!="") < 0) {
  pkc_files <- paste0(pkc_dir, "/", pkc_filenames %>% .[.!=""])
} else {
  pkc_files <- dir(pkc_dir, pattern = pkc_filename_pattern,
                   full.names = TRUE, recursive = TRUE)
}

# Get the sheet name if it's not set.
if(is.null(phenodata_sheet_name) || phenodata_sheet_name=="") {
  # Read in the sample annotation file and get the first sheet name off of it. 
  phenodata_sheet_name <- readxl::excel_sheets(sample_annotation_file)[1] 
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Data import and cleaning ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## ................................................
##
### Data loading ----
##
## ................................................
# Load the data to create a data object using the readNanoStringGeoMxSet function.
# 2024/12/11: We will create a separate object for each PKC module.
data_object_list <- list()
for(pkc_file in pkc_files) {
  # Create the data object (a NanoString GeoMx set) from the input files.
  data_object <- readNanoStringGeoMxSet(dccFiles = dcc_files,
                                        pkcFiles = pkc_file, # Formerly pkc_files, now we only do 1 per object.
                                        phenoDataFile = sample_annotation_file,
                                        phenoDataSheet = phenodata_sheet_name,
                                        phenoDataDccColName = phenodata_dcc_col_name,
                                        protocolDataColNames = protocol_data_col_names,
                                        experimentDataColNames = experiment_data_col_names)
  # Get the name of the current module.
  module <- annotation(data_object) %>% regexPipes::gsub("\\.\\D+$", "")
  
  # Save to the list.
  data_object_list[[module]] <- data_object
  
  # Clean up.
  rm(data_object)
  gc()
}

# For each separate data object: clean and filter the data, and generate neovariables.
for(module in names(data_object_list)) {
  data_object <- data_object_list[[module]]
  
  ## ................................................
  ##
  ### Data cleaning ----
  ##
  ## ................................................
  # Shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
  data_object <- shiftCountsOne(data_object, useDALogic = TRUE)
  
  # Change the column names in the pData to lowercase
  # to match the expected inputs in the NanoString Bioconductor package tools (https://rdrr.io/github/Nanostring-Biostats/GeomxTools/src/R/NanoStringGeoMxSet-qc.R).
  # https://stackoverflow.com/a/51793188/23532435
  # https://stackoverflow.com/questions/69661679/change-multiple-columns-to-lowercase-with-dplyr-difficulty-with-mutate-across-e
  to_lowercase <- c("Area", "Nuclei")
  for(element in to_lowercase) {
    if(element %in% colnames(pData(data_object))) {
      pData(data_object) <- pData(data_object) %>% 
        mutate_at(vars(element), funs(tolower(.)))
    }
  }

  ## ................................................
  ##
  ### Filtering ----
  ##
  ## ................................................
  # Filter observations by the criteria given in filter_vars, if applicable.
  # https://bioconductor.org/packages/release/bioc/vignettes/GeomxTools/inst/doc/Developer_Introduction_to_the_NanoStringGeoMxSet.html
  # filter_vars <- "segment,include,Full ROI,Tumor;Tags,exclude,Stroma"
  
  if(!(is.null(filter_vars) || filter_vars=="")) {
    filter_vars <- filter_vars %>% strsplit(";") %>% unlist
    for(var in filter_vars) {
      var <- var %>% strsplit(",") %>% unlist
      if(!(var[1] %in% colnames(pData(data_object)))) {
        warning(paste0("The variable ", var[1], " is not included in the pData for this data set. Skipping this variable."))
        next
      }
      
      var_values <- var[3:length(var)]
      if((var[2] %>% str_to_lower()) == "include") {
        data_object <- data_object %>% .[,.[[var[1]]] %in% var_values]
      } else if((var[2] %>% str_to_lower()) == "exclude") {
        data_object <- data_object %>% .[,!(.[[var[1]]] %in% var_values)]
      } else {
        warning("Skipping this element of filter_vars. Please enter a value of 'include' or 'exclude' as the second element of filter_vars. And check your spelling!")
        next
      }
      
    }
  }
  
  ## ................................................
  ##
  ### Neovariable generation ----
  ##
  ## ................................................
  # Create the specified neovariables (if provided) and add to the pData.
  if(!is.null(neovariables) && (sum(neovariables != "") > 0)) {
    neovariables <- neovariables %>% .[. != ""]
    for(neovariable in neovariables) {
      # Split into its component parts.
      neovariable_comps <- neovariable %>% strsplit("\\+") %>% unlist
      
      # Create the neovariable.
      neovariable_name <- paste0(neovariable_comps, collapse = "_")
      pData(data_object) <- pData(data_object) %>% tidyr::unite(!!as.name(neovariable_name), neovariable_comps, remove = FALSE, na.rm = FALSE)
    }
  }
  
  ## ................................................
  ##
  ### Save updated object ----
  ##
  ## ................................................
  # Save the updated data object back to the list.
  data_object_list[[module]] <- data_object
  
  # Clean up.
  rm(data_object)
  gc()
  
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Export the NanoStringGeoMxSet object.
saveRDS(data_object_list, paste0(output_dir_rdata, "NanoStringGeoMxSet_raw.rds"))