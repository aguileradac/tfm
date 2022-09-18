#' Load required libraries
#'
#' @description ASWLoadLibraries() loads required libraries.
#'
#' @author David Aguilera Castro
#' 
#' @examples 
#' 
#' ASWLoadLibraries()
#' 

ASWLoadLibraries <- function() {

  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install(version="3.15")
    BiocManager::install('mixOmics', force=T)
  }
  
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load("readxl",
                 "magrittr", 
                 "tibble", 
                 "dplyr", 
                 "SummarizedExperiment",
                 "impute",
                 "PEMM", 
                 "tidyverse",
                 "mixOmics",
                 "class",
                 "e1071",
                 "caret",
                 "doParallel",
                 "MLmetrics", 
                 "mltest",
                 "doBy",
                 "plyr",
                 "heatmaply",
                 "crayon",
                 "bigstatsr",
                 "xtable",
                 "gridExtra"
                 )

}