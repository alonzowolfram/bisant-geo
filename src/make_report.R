#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

# Load the RDS objects with the necessary data (which is just the latest module completed). 
rds_path <- cl_args[5]
latest_module <- readRDS(rds_path)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Make the report -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
message(paste0("Rendering HTML report using the template found at ", rmd_template_file, "."))
# saveRDS(iris, paste0(output_dir_pubs, "test.rds"))
# Set the paths.
# QC - study summary 
pkc_summary_file <- paste0(output_dir_rdata, "pkc_summary_table.rds")
qc_segments_summary_file <- paste0(output_dir_rdata, "qc-segments_summary_table.rds")
ntc_summary_file <- paste0(output_dir_rdata, "ntc_summary_table.rds")
# QC - segments
plot_list_qc_segments_file <- paste0(output_dir_rdata, "qc-segments_plot_list.rds")
# QC - probes
plot_list_qc_probes_file <- paste0(output_dir_rdata, "qc-probes_plot_list.rds")
qc_probes_summary_file <- paste0(output_dir_rdata, "qc-probes_table.rds")
goi_detection_rate_table_file <- paste0(output_dir_rdata, "genes-of-interest_detection-rate_table.rds")
# Normalization
plot_list_normalization_file <- paste0(output_dir_rdata, "normalization_plot_list.rds")
# 16S analysis
plot_list_16s_file <- paste0(output_dir_rdata, "16S-analysis_raw-plots-list.rds")
# Unsupervised analysis
plot_list_unsupervised_clustering_grid_file <- paste0(output_dir_rdata, "plot-list_unsupervised-clustering-grids.rds")
plot_list_cv_heatmap_file <- paste0(output_dir_rdata, "plot-list_cv_heatmaps.rds")
plot_list_16s_score_file <- paste0(output_dir_rdata, "plot-list_16s-score.rds")
# Differential expression analysis
plot_list_diff_exprs_grid_file <- paste0(output_dir_rdata, "LMM-DEG_volcano-plot_grids.rds")
# Pathway analysis
plot_list_pathway_analysis_grid_file <- paste0(output_dir_rdata, "pathway-analysis_grid-arranged-plots-list.rds")
# Immune deconvolution
plot_list_imm_decon_file <- paste0(output_dir_rdata, "immune-deconvolution_plots-list.rds")

# Render.
rmarkdown::render(rmd_template_file,
                  output_file="report.html",
                  output_dir=output_dir_pubs,
                  params=list(qc_segments_summary_file = qc_segments_summary_file,
                              ntc_summary_file = ntc_summary_file,
                              pkc_summary_file = pkc_summary_file,
                              plot_list_qc_segments_file = plot_list_qc_segments_file,
                              plot_list_qc_probes_file = plot_list_qc_probes_file,
                              qc_probes_summary_file = qc_probes_summary_file,
                              goi_detection_rate_table_file = goi_detection_rate_table_file,
                              plot_list_normalization_file = plot_list_normalization_file,
                              plot_list_16s_file = plot_list_16s_file,
                              plot_list_unsupervised_clustering_grid_file = plot_list_unsupervised_clustering_grid_file,
                              plot_list_cv_heatmap_file = plot_list_cv_heatmap_file,
                              plot_list_16s_score_file = plot_list_16s_score_file,
                              plot_list_diff_exprs_grid_file = plot_list_diff_exprs_grid_file,
                              plot_list_pathway_analysis_grid_file = plot_list_pathway_analysis_grid_file,
                              plot_list_imm_decon_file = plot_list_imm_decon_file
                              )
                  )