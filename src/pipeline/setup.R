## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##   License ----
##
#     bisantine-geo is a Snakemake pipeline for processing, running QC on, and analyzing NanoString (Bruker) GeoMx spatial transcriptomics data.
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
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
##                                                                
## R settings ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

options(scipen = 6, digits = 4) # View outputs in non-scientific notation
# memory.limit(30000000)     # This is needed on some PCs to increase memory allowance, but has no impact on Macs

# Prevent regexPipes functions from masking base functions
# https://stackoverflow.com/a/5564654
# It's crazy how many of our problems stem from that lol 8/
grep <- base::grep
grepl <- base::grepl
gsub <- base::gsub

# Load command-line arguments
cl_args <- commandArgs(TRUE)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Functions ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Load required functions
library(tidyverse) # For a data-science focused data "grammar"
## Load helper functions
workflow_system <- cl_args[2]
if(workflow_system=="Nextflow") {
  path <- Sys.getenv("PATH") |> strsplit(":")
  bin_path <- tail(path[[1]], n=1)
  source(file.path(bin_path, "functions/helper_functions.R"))
} else {
  bin_path <- ""
  source("src/functions/helper_functions.R")
}

## Function to add a slash to a directory path if it doesn't have one already
appendSlashToPath <- function(x) {
  ifelse(base::grepl("\\/$", x), x, paste0(x, "/"))
}
## Function to check if variable is NULL or empty
flagVariable <- function(x) {
  return(is.null(x) || sum(x=="") >= length(x))
}

## Function to recursively assign variables
assignVarsEnv <- function(yaml_list, env = .GlobalEnv) { # env = new.env()
  assignRecursive <- function(lst, parent_keys = NULL) {
    for (key in names(lst)) {
      full_key <- paste(c(parent_keys, key), collapse = "_")
      value <- lst[[key]]
      
      if (is.list(value)) {
        assignRecursive(value, c(parent_keys, key))
      } else {
        assign(key, value, envir = env)
      }
    }
  }
  
  assignRecursive(yaml_list)
  # return(env)  # Return the environment so you can use it
}

# Function to validate and process config variables
validateProcessConfig <- function(config_var_metadata) { 
  config_vars <- read.csv(config_var_metadata, stringsAsFactors = FALSE)
  error_msg_list <- list()
  
  for (i in seq_len(nrow(config_vars))) {
    var_name <- config_vars$variable[i]
    required <- as.logical(config_vars$required[i])
    default_value <- config_vars$default_value[i]
    split_string <- as.logical(config_vars$split_string[i])
    
    # Debugging prints
    message(paste("Processing:", var_name))
    message(paste(" - Required:", required))
    message(paste(" - Default Value (raw):", default_value))
    message(paste(" - Exists in .GlobalEnv?", exists(var_name, envir = .GlobalEnv, inherits = FALSE)))
    
    # Convert default_value to the appropriate type
    if (!is.na(default_value) && default_value != "") {
      if (grepl("^\\d+$", default_value)) {  
        # Integer check (only digits, no decimal)
        default_value <- as.integer(default_value)
      } else if (grepl("^\\d+\\.\\d+$", default_value)) {  
        # Decimal check (digits + decimal point)
        default_value <- as.numeric(default_value)
      } else if (default_value %in% c("TRUE", "FALSE")) {
        default_value <- as.logical(default_value)
      }
    }
    
    message(paste(" - Default Value (converted):", default_value, "Type:", typeof(default_value)))
    
    if (exists(var_name, envir = .GlobalEnv, inherits = FALSE)) { 
      value <- get(var_name, envir = .GlobalEnv)
      message(paste(" - Current Value:", value))
      
      if (required && flagVariable(value)) {
        error_msg_list[[var_name]] <- paste("Input missing for", var_name)
      }
      
      if (!required && flagVariable(value) && default_value != "") {
        assign(var_name, default_value, envir = .GlobalEnv)
        message(paste(" - Assigned default:", var_name, "=", get(var_name, envir = .GlobalEnv)))
      }
      
      if (!flagVariable(value) && split_string) {
        assign(var_name, strsplit(value, ",")[[1]] %>% trimws(), envir = .GlobalEnv)
        message(paste(" - Split and assigned:", get(var_name, envir = .GlobalEnv)))
      }
      
    } else {
      if (!required && default_value != "") {
        assign(var_name, default_value, envir = .GlobalEnv)
        message(paste(" - Assigned default:", var_name, "=", get(var_name, envir = .GlobalEnv)))
      }
    }
  }
  
  return(error_msg_list)
}

# Function to generate QC histograms
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

# Function to extract all variables from a formula
extractVariables <- function(formula_str) {
  # Remove "score ~" and split by operators (+ and |)
  formula_vars <- str_remove(formula_str, "\\~") %>%
    str_split("\\+|\\|") %>%
    unlist() %>%
    pipe.gsub("[\\(\\)]", "") %>% 
    str_trim() %>% 
    # Remove "1"s
    pipe.grep("[1]{1}", invert = TRUE, value = TRUE)
  
  # Remove empty elements and return unique variable names
  return(unique(formula_vars[formula_vars != ""]))
}
# Function to extract first fixed effect from a model formula
extractFirstFixedEffect <- function(model_formula) {
  # Convert formula to character and extract terms
  terms <- all.vars(model_formula)
  
  # Find the index of the first random effect (terms with "|")
  random_effect_idx <- which(base::grepl("\\|", attr(terms(model_formula), "term.labels")))
  
  # If there's a random effect, take only fixed effects
  if (length(random_effect_idx) > 0) {
    fixed_effects <- attr(terms(model_formula), "term.labels")[1:(random_effect_idx[1] - 1)]
  } else {
    fixed_effects <- attr(terms(model_formula), "term.labels")
  }
  
  # Return the first fixed effect
  if (length(fixed_effects) > 0) {
    return(fixed_effects[1])
  } else {
    stop("No fixed effects found in the model formula.")
  }
}

# Function to extract all variables from a model formula
extractModelComponents <- function(model_formula) {
  # Extract the terms from the formula
  term_labels <- attr(terms(model_formula), "term.labels")
  
  # Identify random effect terms (contain "|")
  random_terms <- term_labels[base::grepl("\\|", term_labels)]
  
  # Extract fixed effects (everything that is not a random effect)
  fixed_effects <- setdiff(term_labels, random_terms)
  
  # Initialize vectors for random intercepts and slopes
  random_intercepts <- c()
  random_slopes <- c()
  
  # Process random effect terms
  for (rand_term in random_terms) {
    # Split by "|" to separate the grouping factor
    split_term <- strsplit(rand_term, "\\|")[[1]]
    effect_part <- trimws(split_term[1])  # Random effect terms
    grouping_var <- trimws(split_term[2]) # Grouping factor
    
    # Check if it's a random intercept or slope
    if (effect_part == "1") {
      random_intercepts <- c(random_intercepts, grouping_var)
    } else {
      slopes <- unlist(strsplit(effect_part, "\\+"))
      slopes <- trimws(slopes)  # Clean up spaces
      random_slopes <- c(random_slopes, paste(slopes, "by", grouping_var))
    }
  }
  
  # Return as a named list
  return(list(
    fixed_effects = fixed_effects,
    random_intercepts = unique(random_intercepts),
    random_slopes = unique(random_slopes)
  ))
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Parameters from configuration YAML file ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Set the parameters passed from the configuration YAML file
## Read in the config file
library(yaml) # For reading in YAML documents
config <- yaml::read_yaml(cl_args[1])

# Check the required parameters passed from the configuration YAML file based on which module we're running
current_module <- cl_args[3]

## Load the raw variables
config_env <- assignVarsEnv(config)
## Process the variables
config_metadata_path <- "config_variables.csv"
error_msg_list <- validateProcessConfig(config_metadata_path)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Required libraries and functions ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(broom.mixed) # For tidy() function to clean up LMM output
library(cowplot) # For plot_grid
library(emmeans) # Used for post-hoc pairwise comparisons after linear mixed modeling
library(fgsea) # For GSEA
library(GeomxTools) # For NanoString GeoMx stuff
library(ggforce) # `gather_set_data()` function used in study-design QC module
library(ggplotify) # Convert grobs (such as those created by grid) into ggplot objects, so we can display them in Rmd
library(ggpubr) # For annotate_figure(), as_ggplot()
library(ggrepel) # For graphing
library(grid) # For textGrob()
library(gridExtra) # For arrangeGrob / grid.arrange
library(immunedeconv) # One-stop shop for immune deconvolution
library(knitr) # For tables
library(lme4) # Linear mixed models, used in differential expression and differential abundance analyses
library(msigdbr) # Connecting to MSigDB
library(NanoStringNCTools) # For NanoString stuff
library(networkD3) # For Sankey plots
library(officer) # For PowerPoint output
library(openxlsx) # For reading and writing Microsoft Excel files
library(preprocessCore) # For normalize.quantiles()
library(pheatmap) # For heatmaps
library(readxl) # For reading Excel files
library(reshape2) # For melt()
library(Rtsne) # For t-SNE plots
library(scales) # For percents
library(shiny) # For RMD-related stuff. Can't remember exactly
library(SpatialDecon) # For spatial deconvolution
library(stringi) # For string manipulations
library(umap) # For UMAPs
library(vegan) # For ecological diversity indices (Shannon, Simpson, etc.)
# library(devtools)

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
##                                                                
## File/path settings ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
output_dir <- ""
output_dir_config <- ""
output_dir_logs <- ""
output_dir_rdata <- ""
output_dir_tabular <- ""
output_dir_pubs <- ""
output_dir_imgs <- ""

if(workflow_system != "Nextflow") {
  ## Output
  output_dir <- paste0(appendSlashToPath(cl_args[4]))
  ### Create the directory if it doesn't already exist
  if(!dir.exists(output_dir)) dir.create(output_dir)
  ### Create the folder structure within the output_dir
  for(subdir in c("config", "logs", "pubs", "Rdata", "tabular")) {
    subdir_path <- file.path(output_dir, subdir)
    if(!dir.exists(subdir_path)) dir.create(subdir_path)
  }
  
  output_dir_config <- paste0(output_dir, "config/")
  output_dir_logs <- paste0(output_dir, "logs/")
  output_dir_rdata <- paste0(output_dir, "Rdata/")
  output_dir_tabular <- paste0(output_dir, "tabular/")
  output_dir_pubs <- paste0(output_dir, "pubs/")
  
  # Create a subdirectory within `output_dir_pubs` for images
  output_dir_imgs <- paste0(output_dir, "pubs/imgs/")
  if(!dir.exists(output_dir_imgs)) dir.create(output_dir_imgs)
}

rdata_folder <- ifelse(workflow_system=="Nextflow", "", "Rdata/")

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##         
## Parameter check ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# If any of the required parameters are missing, 
# print a message for each one, then stop the pipeline
if(length(error_msg_list) > 0) {
  message("Error: you are missing one or more required parameters. Please see the error messages below.")
  
  for(msg in error_msg_list) {
    message(msg)
  }
  
  stop("Shutting down pipeline due to missing parameters")
}

# There are a couple of parameters we have to handle manually:
# `output_plot_file_types`, `analyte`, `subset_var_16s_levels_manual`, `subset_var_de_levels_manual`, `de_genes_cutoffs`, `imm_decon_methods`, `subset_var_imm_decon_levels_manual`, and `subset_var_tcr_levels_manual`
if(flagVariable(output_plot_file_types)) {output_plot_file_types <- c("eps", "pdf")} else {output_plot_file_types <- output_plot_file_types %>% str_split(",") %>% unlist %>% str_to_lower}
if(flagVariable(analyte) || !(analyte %in% c("RNA", "protein"))) analyte <- "RNA"
if(flagVariable(subset_var_de_levels_manual)) {
  subset_var_de_levels_manual <- rep(list(NA), length(subset_vars_de))
  names(subset_var_de_levels_manual) <- subset_vars_de
} else {
  subset_var_de_levels_manual <- subset_var_de_levels_manual %>% strsplit(";") %>% unlist %>% strsplit(",")
  if(length(subset_var_de_levels_manual)==length(subset_vars_de)) {
    names(subset_var_de_levels_manual) <- subset_vars_de
  } else {
    subset_var_de_levels_manual <- rep(list(NA), length(subset_vars_de))
    names(subset_var_de_levels_manual) <- subset_vars_de
  }
}
if(flagVariable(de_genes_cutoffs)) de_genes_cutoffs <- c(0.25, 0.58)
if(flagVariable(imm_decon_methods)) imm_decon_methods <- c("mcp_counter", "quantiseq")
if(flagVariable(subset_var_imm_decon_levels_manual)) {
  subset_var_imm_decon_levels_manual <- rep(list(NA), length(subset_vars_imm_decon))
  names(subset_var_imm_decon_levels_manual) <- subset_vars_imm_decon
} else {
  subset_var_imm_decon_levels_manual <- subset_var_imm_decon_levels_manual %>% strsplit(";") %>% unlist %>% strsplit(",")
  if(length(subset_var_imm_decon_levels_manual)==length(subset_vars_imm_decon)) {
    names(subset_var_imm_decon_levels_manual) <- subset_vars_imm_decon
  } else {
    subset_var_imm_decon_levels_manual <- rep(list(NA), length(subset_vars_imm_decon))
    names(subset_var_imm_decon_levels_manual) <- subset_vars_imm_decon
  }
}
if(flagVariable(subset_var_16s_levels_manual)) {
  subset_var_16s_levels_manual <- rep(list(NA), length(subset_vars_16s))
  names(subset_var_16s_levels_manual) <- subset_vars_16s
} else {
  subset_var_16s_levels_manual <- subset_var_16s_levels_manual %>% strsplit(";") %>% unlist %>% strsplit(",")
  if(length(subset_var_16s_levels_manual)==length(subset_vars_16s)) {
    names(subset_var_16s_levels_manual) <- subset_vars_16s
  } else {
    subset_var_16s_levels_manual <- rep(list(NA), length(subset_vars_16s))
    names(subset_var_16s_levels_manual) <- subset_vars_16s
  }
}
if(flagVariable(subset_var_tcr_levels_manual)) {
  subset_var_tcr_levels_manual <- rep(list(NA), length(subset_vars_tcr))
  names(subset_var_tcr_levels_manual) <- subset_vars_tcr
} else {
  subset_var_tcr_levels_manual <- subset_var_tcr_levels_manual %>% strsplit(";") %>% unlist %>% strsplit(",")
  if(length(subset_var_tcr_levels_manual)==length(subset_vars_tcr)) {
    names(subset_var_tcr_levels_manual) <- subset_vars_tcr
  } else {
    subset_var_tcr_levels_manual <- rep(list(NA), length(subset_vars_tcr))
    names(subset_var_tcr_levels_manual) <- subset_vars_tcr
  }
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
##                                                                
## Miscellaneous settings ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Normalization names
normalization_names <- c("raw", "Q3-normalized", "background-normalized", "background-subtracted", "background-subtracted + Q3-normalized", "background-subtracted + P90-normalized", "background-subtracted + background-normalized", "quantile-normalized", "housekeeper-normalized")
names(normalization_names) <- c("exprs", "q3_norm", "neg_norm", "bg_sub", "bg_sub_q3", "bg_sub_p90", "bg_sub_neg", "quant", "hk_norm")

# Name of the combined module (WTA+TCR)
combined_module_wta_tcr <- ifelse(!flagVariable(module_tcr), c(main_module, module_tcr) %>% paste(collapse = ","), "") 