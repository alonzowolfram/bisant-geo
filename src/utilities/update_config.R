#!/usr/bin/env Rscript
# https://training.nextflow.io/advanced/structure/#bin

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Setup -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Libraries
library(tidyr)
library(yaml)
library(stringr)

# Read in the variables from the arguments
cl_args <- commandArgs(TRUE)
master_yaml <- cl_args[1] # The main template `config.yaml` file
project_yaml_list <- cl_args[2] %>% str_split(",") %>% unlist # A vector of configuration YAML files to be updated to match `master_yaml`, separated by commas

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
  # Read config files as raw text
  master_lines <- readLines(master_file, warn = FALSE) %>% sapply(trimTrailingWhitespace)
  project_lines <- if (file.exists(project_file)) readLines(project_file, warn = FALSE) %>% sapply(trimTrailingWhitespace) else NULL
  project_lines_non_comment <- project_lines %>% .[!startsWithHash(.)]
  
  # Process project file line by line
  updated_lines <- c()
  for (line in master_lines) {
    # Check if it's a comment or empty/whitespace line
    if(startsWithHash(line) || emptyOrWhitespace(line)) {
      updated_lines <- c(updated_lines, line)
    } else {
      # Check if the variable already exists in the copy (project_file)
      # First extract the key only, not the value
      key <- line %>% regexPipes::gsub(":.+$", "")
      # Now look for this key in project_lines
      # Make sure to check only the non-comment lines
      index <- base::grep(key, project_lines_non_comment)
      
      if(length(index) < 1) {
        # It doesn't exist
        # So add the line from the master file to `updated_lines`
        updated_lines <- c(updated_lines, line)
      } else {
        # It does
        # So add the line from the project file to `updated_lines`
        updated_lines <- c(updated_lines, project_lines_non_comment[index[1]])
      }
      
    }
  }
  
  # Write the updated file
  writeLines(updated_lines, project_file)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Config update -----------------------------------------
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Loop over each project YAML and update to match the template
for (project_yaml in project_yaml_list) {
  mergeConfig(master_yaml, project_yaml)
}