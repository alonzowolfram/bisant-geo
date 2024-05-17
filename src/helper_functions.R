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