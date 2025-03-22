#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Source the setup.R file.
# Note that we cannot use source() directly in Nextflow; see https://stackoverflow.com/questions/76067626/how-to-call-a-function-present-in-a-different-r-script-from-an-r-executable-usin
# and https://community.seqera.io/t/source-another-r-script-in-the-bin-directory/1059
# So we will have to use the workaround suggested above if the workflow system is Nextflow.
cl_args <- commandArgs(TRUE)
workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "setup.R"))
} else {
  bin_path <- ""
  source("src/setup.R") # I guess the path it sources from is the current working directory, not the path the R script lives in.
}

# Automatically list files in each directory for use.
dcc_files <- dir(dcc_dir, pattern = ".dcc$",
                 full.names = TRUE, recursive = TRUE)
pkc_files <- dir(pkc_dir, pattern = pkc_filename_pattern,
                 full.names = TRUE, recursive = TRUE)
# if(!is.null(pkc_filenames) && sum(pkc_filenames!="") < 0) {
#   pkc_files <- paste0(pkc_dir, "/", pkc_filenames %>% .[.!=""])
# } else {
#   pkc_files <- dir(pkc_dir, pattern = pkc_filename_pattern,
#                    full.names = TRUE, recursive = TRUE)
# }

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
### Data loading: all modules ----
##
## ................................................
# Create the data object (a NanoString GeoMx set) from the input files.
# We have to do this because readNanoStringGeoMxSet() complains when we read in the TCR and BiS modules individually.
# We'll split them into their individual modules in the next step.
data_object_all <- readNanoStringGeoMxSet(dccFiles = dcc_files,
                                      pkcFiles = pkc_files,
                                      phenoDataFile = sample_annotation_file,
                                      phenoDataSheet = phenodata_sheet_name,
                                      phenoDataDccColName = phenodata_dcc_col_name,
                                      protocolDataColNames = protocol_data_col_names,
                                      experimentDataColNames = experiment_data_col_names)

## ................................................
##
### Data loading: individual modules ----
##
## ................................................
# For the future: I guess we can split NanoStringGeoMxSet objects by Module pretty easily?
# subset(target_data_object, Module %in% module_tcr)

# Split `data_object_all` into individual modules.
data_object_list <- list()
for(module in data_object_all@annotation %>% regexPipes::gsub("\\.\\D+$", "")) {
  message(paste0("Working on PKC file ", module, ".pkc"))
  
  # Extract the necessary elements.
  assay_data <- data_object_all@assayData$exprs # Matrix or environment containing the DCCs.
  pheno_data <- data_object_all@phenoData # AnnotatedDataFrame containing pData. AnnotatedDataFrame() takes 3 arguments: data = data.frame, varMetadata = data.frame, dimLabels. varMetadata and dimLabels can be missing.
  feature_data <- data_object_all@featureData@data # AnnotatedDataFrame containing fData.
  experiment_data <- data_object_all@experimentData # MIAME (derived from MIAxE; general container for storing experimental metadata; virtual, cannot be instantiated directly)
  annotation <- data_object_all@annotation # Character vector for the PKC file(s)
  dim_labels <- data_object_all@dimLabels # Character vector of length 2: column names to use as labels for the features and samples respectively
  signatures <- data_object_all@signatures # SignatureSet object
  design <- data_object_all@design # [Optional] one-sided formula representing experimental design.
  feature_type <- data_object_all@featureType # Character string indicating if features are "Probe" or "Target" level
  analyte <- data_object_all@analyte # Character string indicating if features are "RNA" or "Protein"
  protocol_data <- data_object_all@protocolData
  # `assay_data`, `annotation` and `feature_data` need to be subset and `feature_data` repackaged into an AnnotatedDataFrame.
  annotation <- paste0(module, ".pkc")
  feature_data <- feature_data %>% dplyr::filter(Module == module)
  assay_data <- assay_data %>% .[feature_data$RTS_ID,,drop=F]
  
  # Create the data object.
  data_object <- NanoStringGeoMxSet(
    assayData = assay_data,
    phenoData = pheno_data,
    featureData = AnnotatedDataFrame(data = feature_data, varMetadata = data_object_all@featureData@varMetadata, dimLabels = data_object_all@featureData@dimLabels, .__classVersion__ = data_object_all@featureData@.__classVersion__),
    experimentData = experiment_data,
    annotation = annotation,
    dimLabels = dim_labels,
    signatures = signatures,
    design = design,
    featureType = feature_type,
    analyte = analyte,
    protocolData = protocol_data,
    check = FALSE # Needed to bypass checks for "valid" NanoStringGeoMxSet objects
  )
  
  # Save to the list.
  data_object_list[[module]] <- data_object
  
  # Clean up.
  rm(data_object)
  gc()
}

## ................................................
##
### Data loading: WTA+TCR ----
##
## ................................................
if(!flagVariable(combined_module)) {
  wta_tcr_pkc_modules <- c(main_module, module_tcr)
  message(paste0("Working on combined WTA+TCR modules ", paste(wta_tcr_pkc_modules, collapse = ", ")))

  # Create the data object (a NanoString GeoMx set) from the input files.
  wta_tcr_pkc_files <- paste0(pkc_dir, "/", wta_tcr_pkc_modules, ".pkc")
  # Check that all of the requisite PKC files exist.
  if(!all(file.exists(wta_tcr_pkc_files))) { warning("One or more of the WTA and/or TCR PKC files is missing. Combined WTA+TCR object will not be created (but separate objects will still be generated.)") } else {
    data_object <- readNanoStringGeoMxSet(dccFiles = dcc_files,
                                          pkcFiles = wta_tcr_pkc_files,
                                          phenoDataFile = sample_annotation_file,
                                          phenoDataSheet = phenodata_sheet_name,
                                          phenoDataDccColName = phenodata_dcc_col_name,
                                          protocolDataColNames = protocol_data_col_names,
                                          experimentDataColNames = experiment_data_col_names)
    
    # Save to the list.
    data_object_list[[combined_module]] <- data_object
    
    # Clean up.
    rm(data_object)
    gc()
  }
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
  to_lowercase <- c("Slide Name", "Scan Name", "Panel", "Roi", "Segment", "Aoi", "Area", "Tags", "Nuclei")
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

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)