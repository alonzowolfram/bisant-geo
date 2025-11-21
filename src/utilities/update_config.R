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
  source("src/pipeline/setup.R") # I guess the path it sources from is the current working directory, not the path the R script lives in.
}

library(yaml)

# Functions
trimTrailingWhitespace <- function(x) {
  base::sub("[ \t]+$", "", x)
}
startsWithHash <- function(x) {
  base::grepl("^[ \t]*#", x)
}
emptyOrWhitespace <- function(x) {
  base::grepl("^\\s*$", x)
  
  # # Example usage:
  # lines <- c("   ", "", "This has content", "\t\t")
  # sapply(lines, emptyOrWhitespace)
}

mergeConfig <- function(master_file, project_file) {
  # Read config files as raw text.
  master_lines <- readLines(master_file, warn = FALSE) %>% sapply(trimTrailingWhitespace)
  project_lines <- if (file.exists(project_file)) readLines(project_file, warn = FALSE) %>% sapply(trimTrailingWhitespace) else NULL
  project_lines_non_comment <- project_lines %>% .[!startsWithHash(.)]
  
  # Process project file line by line.
  updated_lines <- c()
  for (line in master_lines) {
    # Check if it's a comment or empty/whitespace line. 
    if(startsWithHash(line) || emptyOrWhitespace(line)) {
      updated_lines <- c(updated_lines, line)
    } else {
      # Check if the variable already exists in the copy (project_file).
      # First extract the key only, not the value.
      key <- line %>% regexPipes::gsub(":.+$", "")
      # Now look for this key in project_lines.
      # Make sure to check only the non-comment lines. 
      index <- base::grep(key, project_lines_non_comment)
      
      if(length(index) < 1) {
        # It doesn't exist. 
        # So add the line from the master file to `updated_lines`.
        updated_lines <- c(updated_lines, line)
      } else {
        # It does.
        # So add the line from the project file to `updated_lines`.
        updated_lines <- c(updated_lines, project_lines_non_comment[index[1]])
      }
      
    }
  }
  
  # Write the updated file.
  writeLines(updated_lines, project_file)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Config update -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Example usage
master_yaml <- "~/Library/CloudStorage/OneDrive-InsideMDAnderson/PRIME-TR/projects/DSP-P001_bisant-geo/config.yaml"
project_yaml_list <- c("~/Library/CloudStorage/OneDrive-InsideMDAnderson/PRIME-TR/projects/DSP-P001_bisant-geo/profiles/Lon/DSP-044/DSP-044_RNA_config_complete.yaml")

for (project_yaml in project_yaml_list) {
  mergeConfig(master_yaml, project_yaml)
}