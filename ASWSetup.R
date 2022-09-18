#' ASWSetup libraries, options and functions
#'
#' @description Contains global options and functions.
#'
#' @author David Aguilera Castro
#'

# initialization code ####
source("ASWClasses.R")
source("ASWLoadLibraries.R")

ASWLoadLibraries()  

rm('ASWLoadLibraries')

source("ASWPreTreatment.R")
source("ASWFeaturesSelection.R")
source("ASWModels.R")

# Class: ASWSetup ####
setClass("ASWSetup",
         slots = list(
           stepClasses = "list",
           outputFolders = "list",
           defaultOutputFolder = "character",
           defaultOutputFileName = "character",
           datasetFile = "character",
           featuresIDFile = "character",
           methodsAllowed = "list"
         ),
         prototype = list(
           stepClasses = list(
             IMPUTE    = "ASWStepImpute",
             TRANSFORM = "ASWStepTransform",
             DRIFT     = "ASWStepDriftCorrection",
             FILTER    = "ASWStepFilter",
             SCALE     = "ASWStepScale"),
           outputFolders= list(
             IMPUTE    = "output/001_impute/",
             TRANSFORM = "output/002_transform/",
             DRIFT     = "output/003_drift/",
             FILTER    = "output/004_filter/",
             SCALE     = "output/005_scale/"),
           defaultOutputFolder = "output/",
           defaultOutputFileName = "object",
           datasetFile = "dataset/Data_SparklingWines.xlsx",
           featuresIDFile = "dataset/Data_ID.xlsx",
           methodsAllowed = list(
             "impute" = list(KNN="knn",LLS="lls",MF="mf",MEAN="mean",MEDIAN="median"),
             "transform" = list(LOG="log",LOG2="log2",POWER="power"),
             "drift" = list(loess="loess"),
             "filter" = list(filterByRSD="filterByRSD",removeQCs="removeQCs",aggregateSamples="aggregateSamples"),
             "scale" = list(PARETO="pareto", AUTO="auto", RANGE="range", VAST="vast", LEVEL="level")
           )),
)

            
# auxiliar functions ####
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# Colors for console output
primeColor  <- combine_styles(
  make_style("darkgreen"),
  make_style("azure", bg = TRUE)
)

# Format time difference
HMS <- \(x) sprintf('%02d:%02d:%02d', x %/% 60^2, x %% 60^2 %/% 60, floor(x %% 60^2 %% 60))



