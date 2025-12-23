## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Setup ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Source the setup.R file
source("src/pipeline/setup.R")

# Read in the NanoStringGeoMxSet object
target_data_object_list <- readRDS(cl_args[5])
# Set `main_module` if not set already
modules <- names(target_data_object_list)
if(flagVariable(main_module)) main_module <- modules[1]
rm(modules)
# We'll only need the main module for this one
target_data_object <- target_data_object_list[[main_module]]

# Set the normalization method
normalization_method <- normalization_names[names(normalization_names)==normalization_methods[1]]

# Cell profile matrices for SpatialDecon
data("safeTME")
data("safeTME.matches")

# Convert any "NA" in `subset_vars_imm_decon` to "Complete data set"
subset_vars_imm_decon[subset_vars_imm_decon=="NA"] <- "Complete data set"

# Function to perform NNLS for a given y_i = numeric vector of length n, where n = number of protein markers
# (In other words, y_i corresponds to 1 AOI)
do_nnls <- function(y_vec, S) {
  fit <- nnls::nnls(S, y_vec)
  p <- coef(fit)
  # Normalize to sum to 1 (for proportions)
  p / (sum(p) + 1e-6)
}

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Immune deconvolution ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## ................................................
##
### Score calculation ----
##
## ................................................
# Subset the target data object to include only the transcriptome probes
target_data_object_exprs <- subset(target_data_object, Module %in% main_module)
# Extract the expression matrix
exprs_mat <- target_data_object_exprs@assayData[[normalization_methods[1]]]
# If the species is mouse (Mus musculus), we can still run human methods on mouse data, but we just need to convert to their human orthologues
# Also, download the mouse profile matrix using SpatialDecon::download_profile_matrix()
if(species == "Mus musculus") {
  exprs_mat_ortho <- convert_human_mouse_genes(exprs_mat, convert_to = 'human') # convert_human_mouse_genes() in src/helper_functions.R
  system.time(download_profile_matrix(species = "Mouse", age_group = "Adult", matrixname = "ImmuneAtlas_ImmGen")) # https://github.com/Nanostring-Biostats/CellProfileLibrary/blob/master/Mouse/Mouse_datasets_metadata.csv
  # user  system elapsed 
  # 0.138   0.059   1.858 
}

# Loop through imm_decon_methods
imm_decon_res_list <- list()
if(analyte=="protein") {
  
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #### Protein ----
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  for(method in imm_decon_methods) {
    if(!(method %in% c("protein_cell_abundance", "protein_cell_proportions"))) {
      warning(glue::glue("{method} is not currently supported by this pipeline"))
      next
    } else {
      # Error catching: https://stackoverflow.com/a/55937737/23532435
      skip_to_next <- FALSE
      
      # Method: "protein_cell_abundance"
      if(method=="protein_cell_abundance") {
        if(flagVariable(protein_cell_marker_db)) {
          skip_to_next <- TRUE
        } else {
          if(!file.exists(protein_cell_marker_db)) {
            skip_to_next <- TRUE
          }
        }
        
        # If everything is there, make sure it's the correct file format
        # Read in the file and check that it has at least one entry
        message("Checking provided protein marker database")
        {        # https://stackoverflow.com/questions/71082458/why-cant-if-and-else-statements-be-on-separate-lines-in-r
          if(base::grepl("\\.csv$", protein_cell_marker_db)) { imm_marker_db <- read.csv(protein_cell_marker_db, header = T)}
          else if(base::grepl("\\.tsv$", protein_cell_marker_db)) { imm_marker_db <- read.table(protein_cell_marker_db, header = T, sep = "\t") }
          else if(base::grepl("\\.xls.*$", protein_cell_marker_db)) { imm_marker_db <- read_excel(protein_cell_marker_db, col_names = T, na = "NA")}
          else {skip_to_next <- TRUE}
        }
        if(!skip_to_next) if(!(all(c("module", "target") %in% colnames(imm_marker_db)))) skip_to_next <- TRUE
        
        if(skip_to_next) {
          warning(glue::glue("Could not run method {method}. Please provide all required inputs, or check that paths to provided inputs are correct (file exists, file is the correct format)"))
          next
        }
        
        # Perform immune deconvolution
        imm_decon_res <- data.frame()
        for(module in imm_marker_db$module %>% unique) {
          message(glue::glue("Performing immune deconvolution for protein, cell type/module {module}"))
          
          # Get the identity marker proteins
          id_marker_prots <- imm_marker_db %>% dplyr::filter(module == !!module) %>% .[["target"]]
          
          # Calculate the score for each sample
          score_mat <- exprs_mat %>% 
            .[rownames(.) %in% id_marker_prots,,drop=F] %>% 
            # Remove rows with 0 variance (https://stackoverflow.com/questions/50005717/remove-rows-with-zero-variance-in-r)
            .[apply(., 1, var) != 0, ] %>% 
            # 2025/12/19: for now, we won't scale the expression
            # t %>% scale(center = TRUE, scale = TRUE) %>% t %>%
            colSums() %>% 
            as.matrix %>% t
          rownames(score_mat) <- module
          imm_decon_res <- rbind(imm_decon_res, score_mat)
        }
        imm_decon_res <- imm_decon_res %>% as.data.frame %>% tibble::rownames_to_column(var = "cell_type")
        
        if(skip_to_next) {
          warning(glue::glue("An error occurred when trying to run immune deconvolution method {method}. Skipping to the next method"))
          next
        }
        
      } # End "protein_cell_abundance" method
      
      # Method: "protein_cell_proportions"
      if(method=="protein_cell_proportions") {
        if(flagVariable(protein_cell_marker_db)) {
          skip_to_next <- TRUE
        } else {
          if(!file.exists(protein_cell_marker_db)) {
            skip_to_next <- TRUE
          }
        }
        
        # If everything is there, make sure it's the correct file format
        # Read in the file and check that it has at least one entry
        message("Checking provided protein marker database")
        {
          if(base::grepl("\\.csv$", protein_cell_marker_db)) { imm_marker_db <- read.csv(protein_cell_marker_db, header = T)}
          else if(base::grepl("\\.tsv$", protein_cell_marker_db)) { imm_marker_db <- read.table(protein_cell_marker_db, header = T, sep = "\t") }
          else if(base::grepl("\\.xls.*$", protein_cell_marker_db)) { imm_marker_db <- read_excel(protein_cell_marker_db, col_names = T, na = "NA")}
          else {skip_to_next <- TRUE}
        }
        if(!skip_to_next) if(!(all(c("module", "target", "weight") %in% colnames(imm_marker_db)))) skip_to_next <- TRUE
        
        if(skip_to_next) {
          warning(glue::glue("Could not run method {method}. Please provide all required inputs, or check that paths to provided inputs are correct (file exists, file is the correct format)"))
          next
        }
        
        # Create the weights matrix
        targets <- imm_marker_db[["target"]] %>% unique
        modules <- imm_marker_db[["module"]] %>% unique
        weights <- matrix(0, nrow = length(targets), ncol = length(modules), dimnames = list(targets, modules))
        for(module in modules) {
          for(target in targets) {
            ij <- ifelse(nrow(imm_marker_db %>% dplyr::filter(module==!!module & target==!!target)) < 1, 
                         0, 
                         imm_marker_db %>% dplyr::filter(module==!!module & target==!!target) %>% .[["weight"]] %>% .[1] %>% as.numeric)
            weights[target,module] <- ij
          }
        }
        weights[is.na(weights)] <- 0
        
        # Perform immune deconvolution
        # solve `y_i = S*p_i + epsilon` for `p_i`, where
        # `y_i` = vector of marker expression in AOI `i`
        # `S` = weights matrix
        # `p_i` = vector of cell-type abundances in AOI `i` (<- this is what we want)
        exprs_decon <- exprs_mat %>% t %>% .[, rownames(weights), drop = F]
        # Let's double-check the dimensions before performing NNLS
        # where `m` = number of AOIs, `n` = number of protein markers, `p` = number of cell types
        # exprs_decon: `m` × `n`
        # weights: `n` × `p`
        # exprs_decon * weights = `m` × `p` matrix of m AOIs, each AOI with a vector of `p` proportions ✅
        imm_decon_res <- apply(exprs_decon, 1, do_nnls, S = weights)
        rownames(imm_decon_res) <- colnames(weights)
        imm_decon_res <- imm_decon_res %>% as.data.frame %>% tibble::rownames_to_column(var = "cell_type")
        
        if(skip_to_next) {
          warning(glue::glue("An error occurred when trying to run immune deconvolution method {method}. Skipping to the next method"))
          next
        }
        
      } # End "protein_cell_abundance" method
      
      imm_decon_res_list[[method]] <- imm_decon_res
      
    } # End check for valid immune deconvolution method
  } # End for() loop: immune deconvolution methods
  
} else {
  
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #### RNA ----
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  for(method in imm_decon_methods) {
    if(!(method %in% c("quantiseq", "mcp_counter", "xcell", "epic", "abis", "estimate", "spatialdecon", "mmcp_counter"))) {
      warning(glue::glue("{method} is not currently supported by this pipeline. Please note that TIMER, ConsensusTME, and CIBERSORT are not available for this pipeline due to code incompatabilities"))
      next
    } else {
      # Error catching: https://stackoverflow.com/a/55937737/23532435
      skip_to_next <- FALSE
      
      # If the method is mmcp_counter, check if the species is mouse (Mus musculus)
      if(method == "mmcp_counter") {
        if(species != "Mus musculus") {
          warning("mMCP-counter can only be used with mouse data. Skipping to the next method")
          next
        } else {
          # Using a mouse method on mouse data. Do not convert to orthologues
          exprs_mat_effective <- exprs_mat
        }
      } else {
        # If the current method is not mMCP-counter, check the species
        if(species == "Mus musculus") {
          # Using a human method on mouse data. Convert to orthologues
          exprs_mat_effective <- exprs_mat_ortho
        } else {
          # Using a human method on human data. Do not convert to orthologues
          exprs_mat_effective <- exprs_mat
        }
      }
      
      # Special steps needed for SpatialDecon
      if(method == "spatialdecon") {
        # The spatialdecon function takes 3 arguments of expression data:
        #   
        # 1. The normalized data
        # 2. A matrix of expected background for all data points in the normalized data matrix
        # 3. Optionally, either a matrix of per-data-point weights, or the raw data, which is used to derive weights (low counts are less statistically stable, and this allows spatialdecon to down-weight them)
        # https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html
        
        # If the user provided a custom profile matrix, load it
        load_safeTME <- TRUE
        if(!flagVariable(spatial_decon_profile_matrix)) {
          if(file.exists(spatial_decon_profile_matrix)) {
            load_spatial_decon_profile_matrix <- tryCatch(
              {load(spatial_decon_profile_matrix)}, 
              error = function(msg) {return(NA)}
            )
            if(exists("profile_matrix")) load_safeTME <- FALSE
          }
        }
        if(load_safeTME) { 
          message(glue::glue("The default reference matrix, safeTME, will be used with the SpatialDecon method. If this was not your intended reference matrix, you may have either 1) left the `spatial_decon_profile_matrix` field in the configuration YAML file blank or 2) did not fill in a full path to an .RData file containing an object named `profile_matrix`: a matrix with genes as rows and cell types as columns. If you did intend to use safeTME, you can safely ignore this message"))
          profile_matrix <- safeTME 
        }
        
        # Estimate each data point's expected BG from the negative control probes from its corresponding observation
        # Also, the SpatialDecon algorithm will do its own background subtraction
        # so we have to give it data that has been normalized using a non BG-sub method (otherwise we'll get a LinAlg error)
        # we'll go with quant
        negnames <- fData(target_data_object_exprs) %>% 
          dplyr::filter(Negative == TRUE & Module %in% modules_exprs) %>% 
          .$TargetName
        bg <- derive_GeoMx_background(norm = target_data_object_exprs@assayData[["quant"]],
                                      probepool = fData(target_data_object_exprs)$Module,
                                      negnames = negnames)
        
        signif(profile_matrix[seq_len(3), seq_len(3)], 2)
        # heatmap(sweep(safeTME, 1, apply(safeTME, 1, max), "/"),
        #         labRow = NA, margins = c(10, 5))
        
        # Set the cell profile matrix based on the species
        if(species == "Homo sapiens" & !exists("profile_matrix")) {
          cpm <- safeTME
          message("Using safeTME as reference (profile) matrix for SpatialDecon")
        } else {
          cpm <- profile_matrix %>% as.matrix
          message("Using custom reference (profile) matrix for SpatialDecon")
        }
        
        # Run spatial deconvolution
        system.time({
          res <- tryCatch(runspatialdecon(object = target_data_object_exprs,
                                          norm_elt = "quant", # "neg_norm"   "log_norm"   "bg_sub"     "exprs"      "bg_sub_neg" "quant"      "bg_sub_q3"  "q3_norm"   
                                          raw_elt = "exprs",
                                          X = cpm,
                                          align_genes = TRUE#,
                                          #cell_counts = pData(target_data_object_exprs)
          ),
          error = function(e) {skip_to_next <<- TRUE})
        })
        if(class(res) != "logical") imm_decon_res <- res$beta %>% t %>% as.data.frame %>% rownames_to_column("cell_type")
      } else {
        imm_decon_res <- tryCatch(immunedeconv::deconvolute(exprs_mat_effective, method),
                                  error = function(e) {skip_to_next <<- TRUE})
      }
      
      if(skip_to_next) {
        warning(glue::glue("An error occurred when trying to run immune deconvolution method {method}. Skipping to the next method"))
        next
      }
      
      imm_decon_res_list[[method]] <- imm_decon_res
      
    } # End check for valid immune deconvolution method
  } # End for() loop: immune deconvolution methods
} # End else: analyte is RNA

## ................................................
##
### Differential abundance ----
##
## ................................................
# By "differential abundance," we mean _between-sample_ comparisons of immune-cell populations,
# NOT between-cell-type comparisons.
# Per https://omnideconv.org/immunedeconv/articles/immunedeconv.html,
# the following methods allow between-sample comparisons:
# MCP-counter, xCell, TIMER, ConsensusTME, ESTIMATE, ABIS, mMCP-counter, BASE, EPIC, quanTIseq, CIBERSORT abs., seqImmuCC
between_sample_methods <- c("mcp_counter", "xcell", "estimate", "abis", "mmcp_counter", "epic", "quantiseq", "spatialdecon", "protein_cell_abundance")

# Create a list to hold the excluded levels for each grouping variable/first fixed effect
# These will be generated in the loop for the differential analysis
# and then referenced in the visualization loop
excluded_levels_list <- list()

# Check if the formula table has been provided
# If `formula_table_imm_decon_file` is provided, check that it's a valid file
# If not, skip differential analysis
valid_formula_table <- TRUE
if(!flagVariable(formula_table_file_imm_decon)) { 
  if(!file.exists(formula_table_file_imm_decon)) {
    warning("The path to the formula table provided does not exist. Differential analysis will not be performed")
    valid_formula_table <- FALSE
  } else {
    # Read in the file and check that it has at least one entry
    message("Checking provided formula table file")
    if(base::grepl("\\.csv$", formula_table_file_imm_decon)) { formula_table_imm_decon <- read.csv(formula_table_file_imm_decon, header = F)}
    else if(base::grepl("\\.tsv$", formula_table_file_imm_decon)) { formula_table_imm_decon <- read.table(formula_table_file_imm_decon, header = F, sep = "\t") }
    else if(base::grepl("\\.xls.*$", formula_table_file_imm_decon)) { formula_table_imm_decon <- read_excel(formula_table_file_imm_decon, col_names = F, na = "NA")} 
    else {
      warning("The path to the formula table provided does not exist. Differential analysis will not be performed")
      valid_formula_table <- FALSE
    }
  }
} else {
  warning("No formula table file was provided. Differential analysis will not be performed")
  valid_formula_table <- FALSE
}
# Check that the formula table as read in is valid
if(exists("formula_table_imm_decon")) if(is.null(formula_table_imm_decon)) {warning("The path provided does not contain a valid formula table. Differential analysis will not be performed"); valid_formula_table <- FALSE}

# If everything is in place to do differential analysis, run it
if(valid_formula_table) {
  da_res_list <- list()
  
  # If subset variables are not defined, create a dummy subset variable
  subset_vars_imm_decon_orig_flagged <- FALSE # We'll use this in the plotting section, because we need to know if `subset_vars_imm_decon` was originally NULL/empty in the config YAML file
  if(flagVariable(subset_vars_imm_decon)) { subset_vars_imm_decon <- "Complete data set"; subset_vars_imm_decon_orig_flagged <- TRUE } else {subset_vars_imm_decon_orig_flagged <- FALSE}
  
  for(i in 1:nrow(formula_table_imm_decon)) {
    model_number <- i
    da_res_list[[paste0("model_", model_number)]] <- list()
    
    # Extract the formula, whether or not to do all pairwise comparisons, which level to set as the baseline level (if applicable), and any levels of the first fixed effect to exclude
    formula <- formula_table_imm_decon[i,1] %>% as.character
    all_pairwise <- formula_table_imm_decon[i,2]
    baseline_level <- formula_table_imm_decon[i,3] %>% as.character
    excluded_levels <- formula_table_imm_decon[i,4] %>% as.character %>% str_split(";") %>% unlist
    
    # Strip anything before the `~`.
    formula <- formula %>% pipe.gsub("^([[:space:]]|.)*~", "~")
    # Extract variables from formula
    formula_vars <- extractVariables(formula) # Input: string
    # Extract first fixed effect from formula
    first_fixed_effect <- extractFirstFixedEffect(as.formula(formula)) # Input: formula
    # Add the dependent variable
    formula_af <- paste0("score ", formula) %>% as.formula
    
    # Now that we have both the first fixed effect (grouping variable) and associated excluded levels,
    # add to `excluded_levels_list` 
    excluded_levels_list[[first_fixed_effect]] <- excluded_levels
    
    # Ensure Sample is included
    formula_vars <- c("Sample", formula_vars)
    
    # Add the `Sample` column (created from rownames) to pData, then
    # convert all selected columns (except Sample) to factors
    # and subset to include only necessary variables
    # First identify which columns in pData_subset are data frames (e.g. LOQ)
    df_cols <- sapply(pData(target_data_object), is.data.frame)
    
    # Convert only non-data-frame columns (except Sample) to factors
    pData_tmp <- pData(target_data_object) %>%
      tibble::rownames_to_column(var = "Sample") %>% 
      mutate(across(!is.data.frame, ~ if (is.numeric(.)) . else as.factor(.))) 
    # %>% select(all_of(formula_vars))
    # Add "Complete data set" as variable if it doesn't exist already
    if(("Complete data set" %in% subset_vars_imm_decon) & !("Complete data set" %in% colnames(pData_tmp))) pData_tmp$`Complete data set` <- as.factor("Complete data set")
    
    # If `all_pairwise` is FALSE (i.e., we're doing comparisons against a baseline level)
    # then re-order the factor levels of `first_fixed_effect` in `pData_tmp` so that `baseline_level` is the first
    if(is.logical(all_pairwise)) all_pairwise <- T
    if(!all_pairwise & (baseline_level %in% levels(pData_tmp[[first_fixed_effect]]))) {
      non_baseline_levels <- levels(pData_tmp[[first_fixed_effect]]) %>% .[. != baseline_level]
      new_order <- c(baseline_level, non_baseline_levels)
      pData_tmp[[first_fixed_effect]] <- factor(pData_tmp[[first_fixed_effect]], levels = new_order)
    }
    
    # Loop level 2 (deconvolution method)
    for(method in names(imm_decon_res_list)) {
      da_res_list[[paste0("model_", model_number)]][[method]] <- list()
      
      # Check if the method can be used in between-sample comparisons
      if(!(method %in% between_sample_methods)) next
      
      # Convert `immune_data` to long format
      # and merge with `pData_subset`
      immune_long <- imm_decon_res_list[[method]] %>%
        pivot_longer(cols = -cell_type, names_to = "Sample", values_to = "score") %>%
        left_join(pData_tmp, by = "Sample")
      
      # Loop level 3 (subset variable)
      for(subset_var in subset_vars_imm_decon) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete data set"` above
        if(subset_var == first_fixed_effect) next
        
        da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]] <- list()
        
        if(flagVariable(subset_var) | subset_var == "Complete data set") {
          subset_tag <- "All ROIs"
          subset_var <- "Complete data set"
        } else {
          subset_tag <- subset_var
        }
        
        # Get the levels of the current subset_var
        subset_var_levels <- pData_tmp[[subset_var]] %>% levels # Previously as.factor %>% needed because it might be a character vector, but now we convert all non-data-frame columns (except Sample) to factor above
        
        # If subset_var_imm_decon_levels_manual is set, filter subset_var_levels to include only those values
        subset_var_imm_decon_levels_manual_i <- subset_var_imm_decon_levels_manual[[subset_var]]
        if(sum(is.na(subset_var_imm_decon_levels_manual_i)) < length(subset_var_imm_decon_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
          if(sum(subset_var_imm_decon_levels_manual_i %in% subset_var_levels) < 1) {
            warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
          } else { # At least one subset_var_level_manual_i is an actual level of the current subset variable
            subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_imm_decon_levels_manual_i]
          }
        }
        
        # loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
        for(subset_var_level in subset_var_levels) {
          da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]][[subset_var_level]] <- list()
          
          # Get all the samples belonging to the current subset_var_level
          samples <- pData_tmp %>% dplyr::filter(!!as.name(subset_var)==subset_var_level) %>% dplyr::select(Sample) %>% unlist
          # ind <- pData_tmp[[subset_var]] == subset_var_level
          # pData_tmp_sub <- pData_tmp[ind,]
          
          # loop level 5 (cell type) << loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
          for(cell in unique(immune_long$cell_type)) {
            message("")
            message(glue::glue("Calculating differential abundance using formula {formula} | {method} method | subset variable {subset_var} | level {subset_var_level} for cell type {cell} \n"))
            skip_to_next <- FALSE
            
            # Subset data for this cell type, subset variable, and subset variable level
            cell_data <- immune_long %>%
              dplyr::filter((Sample %in% samples) & (cell_type == cell))
            # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
            if(!flagVariable(excluded_levels)) cell_data <- cell_data %>% dplyr::filter(!(!!as.name(first_fixed_effect) %in% excluded_levels))
            
            # Fit the user-defined linear mixed model
            model <- tryCatch(
              expr = lmerTest::lmer(as.formula(formula_af), data = cell_data),
              error = function(e) {
                warning(glue::glue("Error in fitting linear mixed model for model {formula} - {method} - {subset_var} - {subset_var_level} - {cell}: {e$message}"))
                skip_to_next <<- TRUE
              }
            )
            
            if (skip_to_next) next
            message("Successfully calculated differential abundance")
            
            # If `all_pairwise` is TRUE, extract estimated marginal means (EMMs) for pairwise comparisons
            if(all_pairwise) {
              # Create a dynamic formula for emmeans
              emm_formula <- as.formula(paste0("~ ", first_fixed_effect))
              emm <- emmeans::emmeans(model, emm_formula)  # Replace 'Group' with your fixed effect variable
              # Perform all pairwise comparisons (Tukey-adjusted)
              pairwise_results <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey") %>%
                as.data.frame()
              
              # Manually create `model_summary` data frame based on `tidy(model)`
              # Colnames of `model_summary` as created by `tidy(model)`: 
              # [1] "effect"           "fixed_effect"     "baseline"         "term"             "estimate"         "std.error"        "statistic"       
              # [8] "df"               "p.value"          "subset_var"       "subset_var_level" "formula"    
              # Colnames of `pairwise_results`:
              # [1] "contrast" "estimate" "SE"       "df"       "t.ratio"  "p.value"
              # So, to turn `pairwise_results` into `model_summary`:
              # add "effect" = "fixed"
              # add "fixed_effect" = `first_fixed_effect`
              # add "baseline" and "term" by splitting "contrast" by "-" and swapping the order
              # on that note, "estimate" in `pairwise_results` is the inverse of "estimate" in `tidy(model)`, so multiply by -1 to get the correct value
              # "std.error" is just "SE", so just rename "SE"
              # like "estimate", "statistic" should be the inverse of "t.ratio" (rename and multiply by -1)
              # "df" is just "df" (although the values are different between the output of `lmer` and `emmeans`)
              # "p.value" is just "p.value" (although, again, the values are different due to multiple testing correction)
              # "subset_var", "subset_var_level", and "formula" need to be added manually
              model_summary <- pairwise_results %>% 
                dplyr::mutate(effect = "fixed", fixed_effect = first_fixed_effect) %>% 
                tidyr::separate(col = "contrast", into = c("baseline", "term"), sep = " - ") %>% 
                dplyr::mutate(estimate = -1 * estimate, statistic = -1 * t.ratio) %>% 
                dplyr::rename(std.error = SE) %>% 
                dplyr::mutate(subset_var = subset_var, 
                              subset_var_level = subset_var_level, 
                              formula = formula %>% pipe.gsub("~ ", ""),
                              cell_type = cell,
                              method = method)
              model_summary <- model_summary[,c("method", "cell_type", "effect", "fixed_effect", "baseline", "term", "estimate", "std.error", "statistic", "df", "p.value", "subset_var", "subset_var_level", "formula")]
              
            } else { # Comparisons against a single baseline
              
              # Extract fixed effect results and create `model_summary` data frame
              base_value <- rownames(contrasts(pData_tmp[[first_fixed_effect]]))[1] # Get the baseline value
              model_summary <- tryCatch(
                expr = tidy(model, effects = "fixed") %>%
                  dplyr::filter(term != "(Intercept)") %>% # Filter out the intercept
                  dplyr::mutate(metric = metric,
                                baseline = base_value,
                                subset_var = subset_var,
                                subset_var_level = subset_var_level,
                                cell_type = cell, 
                                fixed_effect = first_fixed_effect,
                                method = method, 
                                formula = as.character(formula_af)[3]) %>%  # Add cell type column. baseline = base_value, 
                  dplyr::mutate(term = base::gsub(paste0("^", fixed_effect), "", term)) %>% # Clean up the fixed effect name
                  dplyr::relocate(fixed_effect, .before = term) %>% 
                  dplyr::relocate(baseline, .before = term) %>%
                  dplyr::relocate(cell_type, .before = effect) %>%
                  dplyr::relocate(method, .before = cell_type),
                error = function(e) {
                  warning(glue::glue("Error in extracting fixed effect results for model {formula} - {subset_var} - {subset_var_level}: {e$message}"))
                  skip_to_next <<- TRUE
                }
              ) # End tryCatch() for model summary
              
            } # End else comparisons against a single baseline
            if (skip_to_next) next
            message("Successfully extracted fixed effect results")
            
            # Store results
            da_res_list[[paste0("model_", model_number)]][[method]][[subset_var]][[subset_var_level]][[cell]] <- model_summary
            
          } # End loop level 5 (cell type) << loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
          
        } # End loop level 4 (level of current subset variable) << loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
        
      } # End loop level 3 (subset variable) << loop level 2 (deconvolution method) << loop level 1 (LMM model)
      
      # # Combine results into a data frame
      # results_df <- bind_rows(da_res_list[[paste0("model_", model_number)]][[method]])
      
    } # End loop level 2 (deconvolution methods)
    
  } # End loop level 1 (LMM model)
  
  da_res_df <- bind_rows(rlang::squash(da_res_list)) # `squash` recursively flattens the list.
  
} # End valid formula table check

## ................................................
##
### Visualization ----
##
## ................................................
# https://omnideconv.org/immunedeconv/articles/detailed_example.html
# Rename the samples to reflect their segment and type
pData_tmp <- pData(target_data_object) %>% as.data.frame %>% tibble::rownames_to_column(var = "rownames")
observation_identifiers <- intersect(observation_identifiers, colnames(pData_tmp)) # Make sure the observation identifiers are actually in the pData.
pData_tmp <- pData_tmp %>% tidyr::unite("All ROIs", c(observation_identifiers, rownames), remove = FALSE, na.rm = FALSE, sep = " | ") 
pData_tmp$`All ROIs` <- pData_tmp$`All ROIs` %>% pipe.gsub("\\.dcc", "")
pData_tmp$`Complete data set` <- "Complete data set"
plot_list <- list()

if(subset_vars_imm_decon_orig_flagged) subset_vars_imm_decon <- NULL
if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete data set"

for(method in names(imm_decon_res_list)) {
  plot_list[[method]] <- list()
  df <- imm_decon_res_list[[method]]

  # QuanTIseq, CIBERSORT (absolute), Epic - visualize as stacked bar charts
  # MCP-counter, SpatialDecon - visualize as box/dot/violin plot
  
  for(subset_var in subset_vars_imm_decon) { # You can't loop over a NULL variable, hence the line `if(flagVariable(subset_vars_imm_decon)) subset_vars_imm_decon <- "Complete data set"` above
    plot_list[[method]][[subset_var]] <- list()
    
    if(flagVariable(subset_var)) {
      subset_tag <- "All ROIs"
      subset_var <- "Complete data set"
    } else {
      subset_tag <- subset_var
    }
    
    # Get the levels of the current `subset_var`
    subset_var_levels <- pData_tmp[[subset_var]] %>% as.factor %>% levels # as.factor needed because it might be a character vector.
    
    # If `subset_var_imm_decon_levels_manual` is set, filter `subset_var_levels` to include only those values
    subset_var_imm_decon_levels_manual_i <- subset_var_imm_decon_levels_manual[[subset_var]]
    if(sum(is.na(subset_var_imm_decon_levels_manual_i)) < length(subset_var_imm_decon_levels_manual_i)) { # At least one subset_var_level_manual_i is not NA
      if(sum(subset_var_imm_decon_levels_manual_i %in% subset_var_levels) < 1) {
        warning(glue::glue("None of the manually provided levels for subset variable {subset_var} are present in that variable. All available levels of subset variable {subset_var} will be used"))
      } else { # At least one `subset_var_level_manual_i` is an actual level of the current subset variable
        subset_var_levels <- subset_var_levels %>% .[. %in% subset_var_imm_decon_levels_manual_i]
      }
    }
    
    # loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
    for(subset_var_level in subset_var_levels) {
      plot_list[[method]][[subset_var]][[subset_var_level]] <- list()
      
      # Get all the samples belonging to the current subset_var_level
      ind <- pData_tmp[[subset_var]] == subset_var_level
      pData_tmp_sub <- pData_tmp[ind,]
      
      # Loop level 4 (grouping variable) << loop level 3 (level of current subset variable) << loop level 2 (subset variable) << loop level 1 (model)
      for(grouping_var in grouping_vars_imm_decon) {
        if(grouping_var==subset_var) next
        
        plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]] <- list()
        df2 <- df
        
        # Gather the dataframe for graphing
        if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
          plot_df <- df2 %>% 
            gather(`All ROIs`, fraction, -cell_type)
        } else {
          plot_df <- df2 %>% 
            gather(`All ROIs`, score, -cell_type)
        }
        plot_df$`All ROIs` <- plot_df$`All ROIs` %>% pipe.gsub("\\.dcc", "")
        
        # Add the grouping variable by left_join
        plot_df <- plot_df %>% 
          dplyr::left_join(pData_tmp_sub %>% 
                             dplyr::select(`All ROIs`, !!as.name(grouping_var)), 
                           by = "All ROIs")
        # Replace NAs with character strings.
        plot_df[[grouping_var]] <- plot_df[[grouping_var]] %>% tidyr::replace_na('NA') # https://www.statology.org/replace-na-with-string-in-r/
        # If remove_na_imm_decon is TRUE (or true or True or tRuE or whatever), remove observations with NAs
        if(remove_na_imm_decon | str_to_lower(remove_na_imm_decon)=="true") {
          plot_df <- plot_df %>% dplyr::filter(!!as.name(grouping_var) != "NA")
        }
        # If `excluded_levels` is set, exclude any specified levels from the first fixed effect
        excluded_levels <- excluded_levels_list[[grouping_var]]
        if(!flagVariable(excluded_levels)) plot_df <- plot_df %>% dplyr::filter(!(!!as.name(grouping_var) %in% excluded_levels))
        
        # If the number of unique values of the pData column `grouping_var` > 50, split into multiple groups for graphing
        # https://forum.posit.co/t/diagram-overload-split-data-into-multiple-charts/104355
        unique_values_max <- ifelse(method %in% c("quantiseq", "epic", "cibersort_abs"), 50, 5)
        unique_values <- unique(plot_df[[grouping_var]])
        group <- (1:length(unique_values) %/% unique_values_max) + 1
        names(group) <- unique_values
        plot_df$ChunkingGroup <- group[plot_df[[grouping_var]]]
        
        # Make sure we have enough colors
        n_colors <- plot_df$cell_type %>% unique %>% length
        mycolors <- colorRampPalette(pal_brewer(palette = "Paired")(12))(n_colors) # show_col(pal_brewer()())
        
        # For quantiseq, CIBERSORT absolute (when installed), and EPIC, create a stacked bar plot to show between-cell-type (within-sample) comparisons
        if(method %in% c("quantiseq", "epic", "cibersort_abs")) {
          # Stacked bar chart
          # https://stackoverflow.com/questions/40361800/r-ggplot-stacked-geom-rect
          # ^ for stacked bar charts using geom_rect()
          # See also https://stackoverflow.com/questions/28956442/automatically-resize-bars-in-ggplot-for-uniformity-across-several-graphs-r
          for(group in unique(plot_df$ChunkingGroup)) {
            message(glue::glue("Creating stacked bar graph for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            
            # Create the plot
            plot <- plot_df %>% 
              dplyr::filter(ChunkingGroup==group) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = fraction, fill = cell_type)) + 
              geom_bar(stat = "identity") + scale_fill_manual(values = mycolors) + # , guide = FALSE
              scale_color_manual(values = mycolors) + # , guide = FALSE
              theme_bw() + 
              theme(panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) + 
              scale_x_discrete(limits = rev(levels(imm_decon_res_list[[method]]))) + 
              labs(y = "quantity",
                title = paste0(method, " deconvolution | subset by ", subset_tag, 
                                  " | level ", subset_var_level, 
                                  "\n compartmentalized by ", grouping_var,
                                  " | group ", group))
            if(length(unique_values) > 10) {
              plot <- plot + coord_flip()
            }
            # Add to plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_stacked-bar-graph_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
            }
            
          } # End chunking group for() loop
        } # End if method is one of the stacked barchart ones
        
        # For CIBERSORT (when installed), MCPcounter, xCell, Abis, and Estimate, create a dot plot to show between-sample (within-cell-type) comparisons
        if(method %in% c("cibersort", "mcp_counter", "xcell", "abis", "estimate")) {
          for(group in unique(plot_df$ChunkingGroup)) {
            message(glue::glue("Creating dotplot for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            
            # Create the plot
            plot <- plot_df %>% 
              dplyr::filter(ChunkingGroup==group) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = score, color = cell_type)) +
              geom_point(size = 4) +
              facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
              scale_fill_manual(values = mycolors) + # , guide = FALSE
              scale_color_manual(values = mycolors) + # , guide = FALSE
              coord_flip() +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) +
              labs(title = paste0(method, " deconvolution | subset by ", subset_tag, 
                                  " | level ", subset_var_level, 
                                  "\n compartmentalized by ", grouping_var,
                                  " | group ", group))
            # Add to the plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_dotplot_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 9, units = "in")
            }
            
          } # End chunking group for() loop
        } # End if method is one of the dot-plot ones
        
        # For SpatialDecon, MCPCounter, and protein cell abundance, create a box/violin/dot plot to show between-sample (within-cell-type) comparisons
        if(method %in% c("spatialdecon", "mcp_counter", "protein_cell_abundance")) {
          for(group in unique(plot_df$ChunkingGroup)) {
            message(glue::glue("Creating boxplot for group {group}, grouping variable {grouping_var}, subset variable {subset_var}, level {subset_var_level}, method {method}"))
            
            # Subset to include only the cells in `spatial_decon_cells_plot`, if applicable
            cells_filter <- plot_df$cell_type %>% unique
            if(!flagVariable(spatial_decon_cells_plot) & (sum(spatial_decon_cells_plot %in% (plot_df$cell_type %>% unique)) >= 1)) cells_filter <- intersect(spatial_decon_cells_plot, plot_df$cell_type %>% unique)
            
            # Create the plot
            plot <- plot_df %>% 
              dplyr::filter(ChunkingGroup==group & cell_type %in% cells_filter) %>% 
              ggplot(aes(x = !!as.name(grouping_var), y = score)) +
              geom_boxplot(aes(fill = !!as.name(grouping_var)), 
                           width = 0.3, 
                           lwd = 0.6,
                           alpha = 0.7, 
                           outlier.shape = NA, 
                           staplewidth = 0.3) + 
              geom_jitter(aes(color = !!as.name(grouping_var)), width = 0.15, size = 2, alpha = 0.9) + 
              labs(x = NULL, y = "Abundance score", fill = grouping_var, color = grouping_var) +
              facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
              scale_fill_manual(values = c("#FD6563","#3767A9", "#F2C14E", "#6A4C93")) + # , guide = FALSE
              scale_color_manual(values = c("#FD6563","#3767A9", "#F2C14E", "#6A4C93")) + # , guide = FALSE
              # coord_flip() +
              theme_bw() +
              theme(axis.text.x = element_blank(), # 
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) # +
              # labs(title = paste0(method, " deconvolution | subset by ", subset_tag, 
              #                     " | level ", subset_var_level, 
              #                     "\n compartmentalized by ", grouping_var,
              #                     " | group ", group))
            
            # If `add_boxplot_pvals_imm_decon` is TRUE, 
            # build p-value data frame and add p-values to plot
            # using `stat_pvalue_manual()`
            if(add_boxplot_pvals_imm_decon) {
              # Check if the current `grouping_var` is in `da_res_df` (the data frame with all the differential LMM results)
              if(grouping_var %in% da_res_df$fixed_effect) {
                # Create the p-value data frame
                # p-value data frame: group1, group2, p, y.position, p.label
                pvals_df <- da_res_df %>%
                  dplyr::filter(fixed_effect == grouping_var & subset_var == !!subset_var & subset_var_level == !!subset_var_level) %>%
                  dplyr::select(baseline, term, p.value, cell_type) %>%
                  dplyr::rename(group1 = baseline, group2 = term, p = p.value) %>%
                  # Because the LMM or whatever reformats values with special characters
                  # by enclosing them in parentheses ("()"), we need to strip those parentheses
                  # out of the values, otherwise it will cause weird stuff to happen with the graphing
                  dplyr::mutate(group1 = group1 %>% pipe.gsub("^\\(", "") %>% pipe.gsub("\\)$", ""),
                                group2 = group2 %>% pipe.gsub("^\\(", "") %>% pipe.gsub("\\)$", ""))
                
                # Calculate y.position: slightly above the highest point per facet
                tops <- plot_df %>% # `plot_df` should already have only the samples with the correct `subset_var` and `subset_var_level`
                  group_by(cell_type) %>% # Group by `cell_type`
                  summarise(y.position = max(score, na.rm = TRUE) * 1.05, .groups = "drop") %>%
                  as.data.frame
                
                # Add `y.position` to `pvals_df`
                # The space between brackets should be ~ 10% the range of the points
                ranges <- plot_df %>% # `plot_df` should already have only the samples with the correct `subset_var` and `subset_var_level`
                  group_by(cell_type) %>% 
                  summarise(range = range(score)) %>% 
                  summarise(range = diff(range)) %>% 
                  as.data.frame
                # If any of the ranges are 0, it will throw off things downstream, so replace any 0s with 1
                ranges$range <- ifelse(ranges$range==0,1,ranges$range)
                # Filter to ensure `pvals_df`, `tops`, and `ranges` all have the same cell types
                cell_types_common <- intersect(pvals_df$cell_type %>% unique, plot_df$cell_type %>% unique)
                pvals_df <- pvals_df %>% dplyr::filter(cell_type %in% cell_types_common)
                tops <- tops %>% dplyr::filter(cell_type %in% cell_types_common)
                ranges <- ranges %>% dplyr::filter(cell_type %in% cell_types_common)
                # Calculate bracket spacing
                bracket_spacing <- 0.15 * ranges[,2]; names(bracket_spacing) <- ranges[,1]
                # Because something is going wrong, lemme talk this through lol
                # the highest bracket should be n * bracket_spacing above the lowest bracket
                # where n = (# of pairwise comparisons) - 1
                highest_bracket <- bracket_spacing * (
                  (pvals_df %>% nrow) / length(cell_types_common) - 1
                )
                # Calculate the y-positions of the brackets
                # for(i in 1:length(pvals_df$cell_type %>% unique)) {
                  y.position <- c()
                  for(type in pvals_df$cell_type %>% unique) { # using `pvals_df` because `y.position` will be added to pvals_df, so the cell types should match its order
                    y.position <- c(y.position, (tops %>% dplyr::filter(cell_type==type) %>% .[1,2]) + seq(from = 0, to = highest_bracket %>% .[names(.)==type], by = bracket_spacing %>% .[names(.)==type]))
                  }
                  # vec_length <- length(y.position)
                  # print(glue::glue("{i}: {vec_length}"))
                # }
                pvals_df[["y.position"]] <- y.position
                # Add `label` column 
                # pvals_df[["p.label"]] <- sprintf("p = %.2f", pvals_df$p) # %.2f = 2 decimal places; %.2g = 2 significant figures
                pvals_df <- pvals_df %>% 
                  dplyr::mutate(p.label = dplyr::case_when(
                    p < 0.0001 ~ "****",
                    p < 0.005 ~ "***",
                    p < 0.01 ~ "**",
                    p < 0.05 ~ "*",
                    TRUE ~ sprintf("p = %.2f", p)
                  ))
                
                # Add p-values to plot using `ggpubr::stat_pvalue_manual()`
                plot <- plot + ggpubr::stat_pvalue_manual(
                  pvals_df,
                  y.position = "y.position",
                  label = "p.label",
                  xmin = "group1",
                  xmax = "group2",
                  bracket.size = 0.4,
                  tip.length = 0.01,
                  inherit.aes = FALSE
                ) + 
                  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
                
              } # End check grouping variable is present in LMM results
            } # End check whether user specified to add p-vals to box plots
            
            # Add to the plot list
            plot_list[[method]][[subset_var]][[subset_var_level]][[grouping_var]][[group]] <- plot
            # Save plot to disk
            plot
            file_name <- glue::glue("imm_decon_boxplot_{subset_var}_{subset_var_level}_{method}_{grouping_var}_{group}")
            for(file_type in output_plot_file_types) {
              ggsave(glue::glue("{file_name}.{file_type}"), path = output_dir_imgs, width = 12, height = 16, units = "in")
            }
            
          } # End chunking group for() loop
        } # End if method is SpatialDecon or MCPCounter
        
        # Clean up
        rm(df2)
        gc()
        
      } # End for() loop: grouping variables
      
    } # End for() loop: subset variable levels
    
  } # End for() loop: subset variables
  
} # End for() loop: deconvolution methods

rm(pData_tmp)
gc()

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##                                                                
## Export to disk ----
##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Export deconvolution results as RDS file and Microsoft Excel file
saveRDS(imm_decon_res_list, paste0(output_dir_rdata, "immune-deconvolution_results.rds"))
openxlsx::write.xlsx(imm_decon_res_list, file = paste0(output_dir_tabular, "immune-deconvolution_results_by-method.xlsx"))
# Export differential abundance results as RDS file and Microsoft Excel file
if(exists("da_res_df")) {
  saveRDS(da_res_df, paste0(output_dir_rdata, "immune-deconvolution_differential-abundance_results.rds"))
  openxlsx::write.xlsx(da_res_df, file = paste0(output_dir_tabular, "immune-deconvolution_differential-abundance_results.xlsx"))
}
# Export the raw plots as RDS file
plot_list %>% saveRDS(paste0(output_dir_rdata, "immune-deconvolution_plots-list.rds"))

# Save environment to .Rdata
save.image(paste0(output_dir_rdata, "env_immune_deconvolution.RData"))

# Update latest module completed.
updateLatestModule(output_dir_rdata, current_module)