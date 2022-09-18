#' Main method: Loading setup and datafil
#'
#' @description This file loads ASWSetup and get dataset file.
#'
#' @export
#'
#' @author David Aguilera Castro
#'
##
{
source("ASWSetup.R")

ASW <- new("ASWSetup")

dataset <- loadASWExperiment(ASW@datasetFile)
}

