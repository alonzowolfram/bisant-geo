## Source the setup.R file.
source("src/setup.R")

# Automatically list files in each directory for use.
dcc_files <- dir(dcc_dir, pattern = ".dcc$",
                 full.names = TRUE, recursive = TRUE)
pkc_files <- dir(pkc_dir, pattern = pkc_filename_pattern,
                 full.names = TRUE, recursive = TRUE)

# Get the sheet name if it's not set.
if(is.null(phenodata_sheet_name) || phenodata_sheet_name=="") {
    # Read in the sample annotation file and get the first sheet name off of it. 
    phenodata_sheet_name <- readxl::excel_sheets(sample_annotation_file)[1] 
}

# Load the data to create a data object using the readNanoStringGeoMxSet function.
data_object <- readNanoStringGeoMxSet(dccFiles = dcc_files,
                                      pkcFiles = pkc_files,
                                      phenoDataFile = sample_annotation_file,
                                      phenoDataSheet = phenodata_sheet_name,
                                      phenoDataDccColName = phenodata_dcc_col_name,
                                      protocolDataColNames = protocol_data_col_names,
                                      experimentDataColNames = experiment_data_col_names)

# Shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
data_object <- shiftCountsOne(data_object, useDALogic = TRUE)

# Export the NanoStringGeoMxSet object.
saveRDS(data_object, paste0(output_dir_rdata, "NanoStringGeoMxSet_raw.rds"))

## Initialize PowerPoint file.
if(!is.null(ppt_template_file) && ppt_template_file != "") {
    pptx <- read_pptx(ppt_template_file)
} else {
    pptx <- read_pptx()
}

# num_template_slides <- length(pptx) # We'll use this later when we move all the template slides to the back.
# Add title slide.
the_date <- Sys.Date() %>%
  format('%B %d, %Y') %>%
  as.character()
# https://www.stat.berkeley.edu/~s133/dates.html
pptx <- pptx %>%
  officer::add_slide(layout = "Title Slide", master = "Office Theme") %>%
  officer::ph_with(value = paste0(project_name, " GeoMx data analysis report"), 
                   location = ph_location_label(ph_label = "Title 1")) %>% # For future reference, see also https://stackoverflow.com/questions/58508859/r-officer-package-how-to-specify-a-certain-placeholder-when-there-are-multiple
  officer::ph_with(value = the_date,
                   location = ph_location_label(ph_label = "Subtitle 2"))
print(pptx, file.path(output_dir_pubs, ppt_output_file))

# # Get the layout summary and properties of the template.
# layout_summary(pptx)
# layout_properties(pptx)

