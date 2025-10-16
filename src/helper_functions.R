#' Export RDS file containing name of latest completed module:
#' 
#' @param output_path Directory to which Rds file should be exported
#' @param module The name of the current module that was just completed
#' @return A Boolean
#' @export
updateLatestModule <- function(output_path, module) {
  file_name <- "latest_rule.Rds"
  latest_rule <- module
  saveRDS(latest_rule, paste0(output_path, file_name))
  
  return(TRUE)
}

#' Calculate coefficient of variation.
#' 
#' @param x A vector.
#' @examples
#' calc_CV(vector)
calc_CV <- function(x) {sd(x) / mean(x)}

#' Save a plot to EPS.
#' 
#' @param plot A plot object (ggplot, etc.)
#' @param path A character string.
#' @param width A numeric value.
#' @param height A numeric value.
#' @returns A boolean.
#' @examples
#' saveEPS(plot, "path/to/plot.eps", width = 12, height = 12)
saveEPS <- function(plot, path, width, height) {
  setEPS()
  postscript(path, width = width, height = height)
  plot
  dev.off()
  
  return(TRUE)
}

#' Save a plot to PNG.
#' 
#' @param plot A plot object (ggplot, etc.)
#' @param path A character string.
#' @param width A numeric value.
#' @param height A numeric value.
#' @returns A boolean.
#' @examples
#' savePNG(plot, "path/to/plot.png", width = 12, height = 12, units = "in", res = 300)
savePNG <- function(plot, path, width, height, units, res) {
  png(path, width = width, height = height, units = units, res = res)
  # https://stackoverflow.com/questions/9206110/using-png-function-not-working-when-called-within-a-function
  # ^Why we need to call print().
  print({plot})
  dev.off()
  
  return(TRUE)
}

# https://github.com/omnideconv/immunedeconv/blob/HEAD/R/mouse_deconvolution_methods.R
#' This function converts the mouse gene symbols into corresponding human ones, and vice versa.
#'
#' This function relies on the `biomaRt`` package and connects to the ENSEMBL repository
#'  to retrieve the gene symbols. If ENSEMBL cannot be reached, another solution will be
#'  used. Since it is memory intensive, users can choose not to run it.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples.
#'    Gene symbols must be the rownames of the matrix. If genes map to multiple gene symbols, the median
#'    expression will be returned.
#'    Can also be a m x 1 vector of the gene names. IN this case, only the converted genes will
#'    be returned.
#' @param mirror the ensembl mirror to use. Possible choices are 'www' (default),
#'    'uswest', 'useast', 'asia'
#' @param other_annot boolean, wether to run the alternative conversion method
#' @param convert_to one of 'human' or 'mouse'. Specifies the organism of the orthologs to look for
#' @return the same matrix, with the counts for the corresponding human genes.
#'    This matrix can directly be used with the immunedeconv methods. A message
#'    will display the ratio of original genes which were converted.
#'
#' @export
convert_human_mouse_genes <- function(gene_expression_matrix, mirror = "www",
                                      other_annot = TRUE, convert_to = c("human", "mouse")) {
  vect <- FALSE
  if (!is.vector(gene_expression_matrix)) {
    gene.names <- rownames(gene_expression_matrix)
    gene_expression_matrix <- as.data.frame(gene_expression_matrix)
    gene_expression_matrix$gene_name <- gene.names
  } else {
    gene_expression_matrix <- data.frame(
      "gene_name" = gene_expression_matrix,
      "X" = rep(0, length(gene_expression_matrix))
    )
    vect <- TRUE
  }
  
  genes.retrieved <- NULL
  tryCatch(
    expr = {
      human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror)
      mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = mirror)
      
      if (convert_to == "human") {
        mart.use <- mouse
        mart.link <- human
        attr <- "mgi_symbol"
        attr.link <- "hgnc_symbol"
      } else {
        mart.use <- human
        mart.link <- mouse
        attr <- "hgnc_symbol"
        attr.link <- "mgi_symbol"
      }
      
      genes.retrieved <<- getLDS(
        attributes = c(attr),
        filters = attr, values = gene.names,
        mart = mart.use, attributesL = c(attr.link), martL = mart.link, uniqueRows = T
      )
      
      if (convert_to == "human") {
        newGenes.counts <<- gene_expression_matrix %>%
          left_join(., genes.retrieved, by = c("gene_name" = "MGI.symbol")) %>%
          select(., -c("gene_name")) %>%
          select(., c("HGNC.symbol", everything()))
      } else {
        newGenes.counts <<- gene_expression_matrix %>%
          left_join(., genes.retrieved, by = c("gene_name" = "HGNC.symbol")) %>%
          select(., -c("gene_name")) %>%
          select(., c("MGI.symbol", everything()))
      }
    },
    error = function(e) {
      if (other_annot) {
        print("Cannot connect to ENSEMBL. Using alternative method. This will take some time.")
        # Code adapted from: https://support.bioconductor.org/p/129636/#9144606
        
        mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
        
        find_corr_gene <- function(gene, mouse_human_genes_df, convert_to = c("human", "mouse")) {
          if (convert_to == "human") {
            orgn.name <- "mouse, laboratory"
            new.orgn <- "human"
          } else {
            orgn.name <- "human"
            new.orgn <- "mouse, laboratory"
          }
          
          class_key <- (mouse_human_genes_df %>%
                          filter(Symbol == gene & Common.Organism.Name == orgn.name))[["DB.Class.Key"]]
          
          if (!identical(class_key, integer(0))) {
            output <- NULL
            new_genes <- (mouse_human_genes_df %>% filter(DB.Class.Key == class_key & Common.Organism.Name == new.orgn))[, "Symbol"]
            
            for (new_gene in new_genes) {
              output <- append(output, new_gene)
            }
            
            if (!is.null(output)) {
              return(data.frame(
                "new_gene" = output,
                "old_gene" = gene
              ))
            }
          }
        }
        
        genes.retrieved <- map_dfr(gene.names, function(x) find_corr_gene(x, mouse_human_genes, convert_to))
        
        newGenes.counts <<- gene_expression_matrix %>%
          left_join(., genes.retrieved, by = c("gene_name" = "old_gene")) %>%
          select(., -c("gene_name")) %>%
          select(., c("new_gene", everything()))
      }
    }
  )
  
  colnames(newGenes.counts)[1] <- "gene_name"
  newGenes.counts <- newGenes.counts[!(is.na(newGenes.counts$gene_name)), ] %>%
    group_by(gene_name) %>%
    summarise_all(median)
  
  newGenes.counts <- as.data.frame(newGenes.counts)
  
  rownames(newGenes.counts) <- newGenes.counts$gene_name
  newGenes.counts <- select(newGenes.counts, -c("gene_name"))
  
  fraction <- 100 * (nrow(newGenes.counts) / nrow(gene_expression_matrix)) %>%
    round(., 1)
  
  message(paste0("ATTENTION: Only the ", fraction, "% of genes was maintained"))
  
  if (!vect) {
    return(newGenes.counts)
  } else {
    return(rownames(newGenes.counts))
  }
}

#' Function to truncate strings
#'
#' @param vec A vector of strings
#' @param max_length Maximum length a string can be before getting truncated
#' @param trunc_length Length to truncate a string to that exceeds max_length
#' @return vector with truncated elements
#'
#' @export
truncate_strings <- function(vec, max_length = 35, trunc_length = 31) {
  sapply(vec, function(x) {
    if (nchar(x) > max_length) {
      paste0(substr(x, 1, trunc_length), " ...")
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

#' Function to check if a vector of RGB colors are valid
#' 
#' @param x A vector of strings
#' @return A vector of Boolean values
#' 
#' @export
areColors <- function(x) {
  vapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  }, FUN.VALUE = TRUE) # https://stackoverflow.com/questions/51657605/fun-value-argument-in-vapply
}