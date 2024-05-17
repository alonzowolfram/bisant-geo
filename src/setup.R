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
normalization_methods <- experiment$normalization$normalization_methods
ann_of_interest <- experiment$normalization$ann_of_interest
### Unsupervised analysis
compartment_vars <- experiment$unsupervised$compartment_vars
heatmap_ann_vars <- experiment$unsupervised$heatmap_ann_vars
### Linear mixed models/differential expression
lmm <- experiment$lmm
random_slope <- lmm$random_slope
test_vars <- lmm$test_vars
random_intercept_vars <- lmm$random_intercept_vars
# random_slope_vars <- lmm$random_slope_vars
subset_vars <- lmm$subset_vars
cv_cutoff <- lmm$cv_cutoff
n_top_genes <- lmm$n_top_genes
de_genes_cutoffs <- lmm$de_genes_cutoffs
### Pathway analysis
pathway_analysis <- experiment$pathway_analysis
species <- pathway_analysis$species
pathway_table_file <- pathway_analysis$pathway_table_file
individual_pathways <- pathway_analysis$individual_pathways
n_max_pathways <- pathway_analysis$n_max_pathways
### Immune deconvolution
immune_deconvolution <- experiment$immune_deconvolution
path_to_lm22 <- immune_deconvolution$path_to_lm22
path_to_cibersort <- immune_deconvolution$path_to_cibersort
imm_decon_methods <- immune_deconvolution$imm_decon_methods
individual_identifier <- immune_deconvolution$individual_identifier
compartment_identifier <- immune_deconvolution$compartment_identifier
segment_identifier <- immune_deconvolution$segment_identifier

### Miscellaneous
random_seed <- experiment$misc$random_seed

## ---------------------------

# Check the required parameters passed from the configuration YAML file based on which module we're running.
current_module <- cl_args[3]
## Data
if(current_module=="data_import_cleaning") {
  if(is.null(config$project$technical$path_to_regexPipes)) stop("Please provide a path to the regexPipes package in the configuration YAML file (path_to_regexPipes).")
  if(is.null(config$data$dcc_dir)) stop("In the configuration YAML file, please provide the absolute path to the directory containing your DCC files.")
  if(is.null(config$data$pkc_dir)) stop("In the configuration YAML file, please provide the absolute path to the directory containing your PKC files.")
  if(is.null(config$data$sample_annotation_file)) stop("In the configuration YAML file, please provide the absolute path to your sample annotation Excel file.")
  if(list(NULL) %in% config$experiment$annotation) stop("In the configuration YAML file, please provide values for all experimental annotation column name settings.")
}
## Experiment
if(current_module=="qc_segments") if(list(NULL) %in% config$experiment$segment_qc) stop("In the configuration YAML file, please provide values for all segment QC settings to run the segment QC module.")
if(current_module=="qc_probes") if(list(NULL) %in% config$experiment$probe_qc) stop("In the configuration YAML file, please provide values for all probe QC settings to run the probe QC module.")
if(current_module=="differential_expression_analysis") if(list(NULL) %in% list(config$experiment$lmm$random_slope, config$experiment$lmm$test_vars, config$experiment$lmm$random_intercept_vars)) stop("In the configuration YAML file, please provide values for all linear mixed model settings (random_slope, test_vars, random_intercept_vars) to run the differential expression module.")
if(current_module=="pathway_analysis") if(is.null(species) | species=="") stop("In the configuration YAML file, please provide a value for the 'species' variable to run the pathway analysis module.")
if(current_module=="immune_deconvolution") if(list(NULL) %in% list(path_to_lm22, path_to_cibersort)) stop("In the configuration YAML file, please provide values for the LM22 path and the CIBERSORT.R path.")
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
if(is.null(normalization_methods) || normalization_methods == "") {
  normalization_methods <- c("quantile")
} else {
  normalization_methods <- normalization_methods %>% strsplit(",") %>% unlist
}
## Unsupervised analysis
if(is.null(compartment_vars) || compartment_vars == "") {
  compartment_vars <- ann_of_interest
} else {
  compartment_vars <- compartment_vars %>% strsplit(",") %>% unlist
}
if(!is.null(heatmap_ann_vars) & heatmap_ann_vars != "") {
  heatmap_ann_vars <- heatmap_ann_vars %>% strsplit(",") %>% unlist
}
## Differential expression analysis.
if(is.null(subset_vars) || subset_vars == "") {
  subset_vars <- NA
} else {
  subset_vars <- subset_vars %>% strsplit(",") %>% unlist
}
if(is.null(n_top_genes) || n_top_genes == "" || !is.integer(n_top_genes)) {
  n_top_genes <- 15
}
if(is.null(de_genes_cutoffs) || de_genes_cutoffs == "") {
  de_genes_cutoffs <- c(0.25, 0.58)
} else {
  # Check if there is a comma.
  if(base::grepl(",", de_genes_cutoffs)) {
    # There is a comma.
    # Split and check.
    de_genes_cutoffs <- de_genes_cutoffs %>% strsplit(",") %>% unlist %>% as.numeric() %>% .[1:2] # Keep only the first two elements, in case the user entered more than 1 comma. 
    for(i in length(de_genes_cutoffs)) {
      if(i == 1) {
        if(is.na(de_genes_cutoffs[i])) de_genes_cutoffs[i] <- 0.25 # User did not provide FDR, did provide LFC.
      } else {
        if(is.na(de_genes_cutoffs[i])) de_genes_cutoffs[i] <- 0.58 # User did not provide LFC, did provide FDR. 
      }
    }
  } else {
    # No comma <= user provided only the FDR cutoff.
    de_genes_cutoffs <- c(de_genes_cutoffs, 0.58)
  }
}
## Immune deconvolution.
if(is.null(imm_decon_methods) || imm_decon_methods == "") {
  imm_decon_methods <- c("cibersort", "cibersort_abs", "mcp_counter")
} else {
  imm_decon_methods <- imm_decon_methods %>% strsplit(",") %>% unlist
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
# random_slope_vars <- random_slope_vars %>% as.character %>% strsplit(",") %>% unlist
individual_pathways <- individual_pathways %>% as.character %>% strsplit(",") %>% unlist

# Create full names for the normalization methods.
normalization_names <- c("raw", "Q3-normalized", "background-normalized", "background-subtracted", "background-subtracted + Q3-normalized", "background-subtracted + background-normalized", "quantile-normalized")
names(normalization_names) <- c("exprs", "q3_norm", "neg_norm", "bg_sub", "bg_sub_q3", "bg_sub_neg", "quant")

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
library(ggpubr) # For annotate_figure(), as_ggplot().
library(grid) # For textGrob().
library(gridExtra) # Not sure, but I'm using so many plot-related packages, why not just throw in another one.
library(umap) # For UMAPs.
library(Rtsne) # For t-SNE plots.
library(ggrepel) # For graphing. 
library(fgsea) # For GSEA.
library(msigdbr) # Connecting to MSigDB.
library(immunedeconv) # One-stop shop for immune deconvolution.
library(e1071) # Required by CIBERSORT.
library(parallel) # Required by CIBERSORT.
library(preprocessCore) # Required by CIBERSORT.
library(GenomicRanges) # https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
library(rtracklayer) # https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
library(Rsamtools) # https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
library(pheatmap) # For heatmaps.
library(openxlsx) # For reading and writing Microsoft Excel files.
# library(caret) # For nearZeroVar().
install.packages(path_to_regexPipes, repos = NULL, type = "source")
library(regexPipes) # For pipe-friendly version of base R's regex functions.
# library(devtools)
# devtools::install_github("jdstorey/qvalue") # Dependency for glmmSeq. https://github.com/StoreyLab/qvalue
# devtools::install_github("myles-lewis/glmmSeq") 
# library(glmmSeq) # General linear mixed models. 

## Load helper functions.
source("src/helper_functions.R")
