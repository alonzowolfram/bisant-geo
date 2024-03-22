## ---------------------------
#     bisantine-geo is a Snakemake pipeline for processing, running QC on, and analyzing NanoString GeoMx spatial transcriptomics data.
#     Copyright (C) 2024 Lon Fong.

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
## ---------------------------

# Set R settings.

options(scipen = 6, digits = 4) # View outputs in non-scientific notation.
# memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on Macs.

## ---------------------------

# Load required functions.
library(tidyverse) # For a data-science focused data "grammar".
## Function to add a slash to a directory path if it doesn't have one already.
appendSlashToPath <- function(x) {
  ifelse(base::grepl("\\/$", x), x, paste0(x, "/"))
}

## Function to generate QC histograms.
QCHistogram <- function(assay_data = NULL,
                        annotation = NULL,
                        fill_by = NULL,
                        thr = NULL,
                        scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

## ---------------------------

# Set the parameters passed from the configuration YAML file.
## Read in the config file. 
cl_args <- commandArgs(TRUE)
library(yaml) # For reading in YAML documents.
config <- yaml::read_yaml(cl_args[1])
output_dir <- cl_args[2]

## Load the raw variables.
## Project
### Meta
meta <- config$project$meta
project_name <- meta$project_name
run_name <- meta$run_name
### Technical
path_to_regexPipes <- config$project$technical$path_to_regexPipes
## Data
data <- config$data
dcc_dir <- data$dcc_dir
pkc_dir <- data$pkc_dir
pkc_filename_pattern <- data$pkc_filename_pattern
pkc_filenames <- data$pkc_filenames
sample_annotation_file <- data$sample_annotation_file
phenodata_sheet_name <- data$phenodata_sheet_name
ppt_template_file <- data$ppt_template_file
previous_run_out_dir <- data$previous_run_out_dir
## Outputs
# (output_dir read in from command line arguments)
## Experiment
experiment <- config$experiment
### Annotation
annotation <- experiment$annotation
phenodata_dcc_col_name <- annotation$phenodata_dcc_col_name
protocol_data_col_names <- annotation$protocol_data_col_names
experiment_data_col_names <- annotation$experiment_data_col_names
### Segment QC
segment_qc <- experiment$segment_qc
min_segment_reads <- segment_qc$min_segment_reads
percent_trimmed <- segment_qc$percent_trimmed
percent_stitched <- segment_qc$percent_stitched
percent_aligned <- segment_qc$percent_aligned
percent_saturation <- segment_qc$percent_saturation
min_negative_count <- segment_qc$min_negative_count
max_ntc_count <- segment_qc$max_ntc_count
min_nuclei <- segment_qc$min_nuclei
min_area <- segment_qc$min_area
### Probe QC
probe_qc <- experiment$probe_qc
min_probe_ratio <- probe_qc$min_probe_ratio
percent_fail_grubbs <- probe_qc$percent_fail_grubbs
loq_sd_cutoff <- probe_qc$loq_sd_cutoff
min_loq <- probe_qc$min_loq
gene_detection_rate <- probe_qc$gene_detection_rate
percent_of_segments <- probe_qc$percent_of_segments
### Normalization
normalization_method <- experiment$normalization$normalization_method
ann_of_interest <- experiment$normalization$ann_of_interest
### Unsupervised analysis
compartment_vars <- experiment$unsupervised$compartment_vars
### Linear mixed models/differential expression
lmm <- experiment$lmm
random_slope <- lmm$random_slope
test_vars <- lmm$test_vars
random_intercept_vars <- lmm$random_intercept_vars
random_slope_vars <- lmm$random_slope_vars
### Miscellaneous
random_seed <- experiment$misc$random_seed

# Check the required parameters passed from the configuration YAML file based on which module we're running.
current_module <- cl_args[3]
## Data
if(current_module=="data_import_cleaning") {
  if(is.null(config$project$technical$path_to_regexPipes)) stop("Please provide a path to the regexPipes package in the config.yaml file (path_to_regexPipes).")
  if(is.null(config$data$dcc_dir)) stop("Please provide the absolute path to the directory containing your DCC files.")
  if(is.null(config$data$pkc_dir)) stop("Please provide the absolute path to the directory containing your PKC files.")
  if(is.null(config$data$sample_annotation_file)) stop("Please provide the absolute path to your sample annotation Excel file.")
  if(list(NULL) %in% config$experiment$annotation) stop("Please provide values for all experimental annotation column name settings.")
} 
## Experiment
if(current_module=="qc_segments") if(list(NULL) %in% config$experiment$segment_qc) stop("Please provide values for all segment QC settings to run the segment QC module.")
if(current_module=="qc_probes") if(list(NULL) %in% config$experiment$probe_qc) stop("Please provide values for all probe QC settings to run the probe QC module.")
if(current_module=="differential_expression_analysis") if(list(NULL) %in% config$experiment$lmm) stop("Please provide values for all linear mixed model settings to run the differential expression module.")
if(is.null(ann_of_interest) | ann_of_interest=="") stop("Please provide values for ann_of_interest.")

# Set values for optional/program-set parameters from the configuration YAML file.
## Data
if(!is.null(previous_run_out_dir) && previous_run_out_dir != "") {
  # previous_run_out_dir provided, check if there's a PowerPoint.
  pptx_files <- list.files(file.path(previous_run_out_dir, "pubs"), full.names = TRUE) %>% regexPipes::grep("\\.pptx", value = TRUE)
  if(length(pptx_files) > 0) {
    # We do have a PowerPoint file, use the first one available. 
    ppt_template_file <- pptx_files[1]
  }
}
if(is.null(pkc_filename_pattern) || pkc_filename_pattern=="") {
  pkc_filename_pattern==".pkc$"
}
if(!is.null(pkc_filenames) && pkc_filenames != "") {
  pkc_filenames <- pkc_filenames %>% strsplit(",") %>% unlist
}
## Output
### Create the directory if it doesn't already exist. 
if(!dir.exists(output_dir)) dir.create(output_dir)
### Create the folder structure within the output_dir.
for(subdir in c("config", "logs", "pubs", "Rdata", "tabular")) {
  subdir_path <- file.path(output_dir, subdir)
  if(!dir.exists(subdir_path)) dir.create(subdir_path)
}
output_dir_config <- paste0(output_dir, "config/")
output_dir_logs <- paste0(output_dir, "logs/")
output_dir_rdata <- paste0(output_dir, "Rdata/")
output_dir_tabular <- paste0(output_dir, "tabular/")
output_dir_pubs <- paste0(output_dir, "pubs/")

# Experiment
## Normalization
normalization_method <- ifelse((is.null(normalization_method) || normalization_method == ""), "q_norm", normalization_method)
## Unsupervised analysis
if(is.null(compartment_vars) || compartment_vars == "") {
  compartment_vars <- ann_of_interest
} else {
  compartment_vars <- compartment_vars %>% strsplit(",") %>% unlist
}
## Miscellaneous
random_seed <- ifelse((is.null(random_seed) || random_seed==""), 1026, random_seed) # E.g. 1026. 

# Clean and format parameter values.
## Input file and folder paths
dcc_dir <- dcc_dir %>% appendSlashToPath() # This is the directory in which the DCC files are located.
pkc_dir <- pkc_dir %>% appendSlashToPath() # This is the directory in which the PKC files are located. 
protocol_data_col_names <- protocol_data_col_names %>% strsplit(",") %>% unlist # E.g. "Aoi,Roi" -> c("Aoi", "Roi")
experiment_data_col_names <- experiment_data_col_names %>% strsplit(",") %>% unlist # E.g. "Panel" -> c("Panel")

## Output file and folder paths
# project_name_string <- ifelse((is.null(project_name) || project_name==""), "", paste0(project_name, "_"))
# run_name_string <- ifelse((is.null(run_name) || run_name==""), "", paste0(run_name, "_"))
ppt_output_file <- paste0("GeoMx-analysis_PowerPoint-report.pptx")
rds_output_file <- list()
for(stage in c("raw", "segment_qc", "probe_qc", "normalized", "unsupervised_clustering")) {
  rds_output_file[[stage]] <- paste0("NanoStringGeoMxSet", "_", stage, ".rds")
}

## Experimental
random_slope <- random_slope %>% as.character %>% strsplit(",") %>% unlist
test_vars <- test_vars %>% as.character %>% strsplit(",") %>% unlist
random_intercept_vars <- random_intercept_vars %>% as.character %>% strsplit(",") %>% unlist
random_slope_vars <- random_slope_vars %>% as.character %>% strsplit(",") %>% unlist

## ---------------------------

## Load required libraries.
library(shiny) # For RMD-related stuff. Can't remember exactly. 
library(knitr) # For tables.
library(officer) # For PowerPoint output.
library(readxl) # For reading Excel files.
library(NanoStringNCTools) # For NanoString stuff.
library(GeomxTools) # For NanoString GeoMx stuff. 
library(stringi) # For string manipulation.
library(ggforce) # I have no idea.
library(scales) # For percents.
library(reshape2) # For melt.
library(cowplot) # For plot_grid.
library(umap) # For UMAPs.
library(Rtsne) # For t-SNE plots.
library(ggrepel) # For graphing. 
install.packages(path_to_regexPipes, repos = NULL, type = "source")
library(regexPipes) # For pipe-friendly version of base R's regex functions.
