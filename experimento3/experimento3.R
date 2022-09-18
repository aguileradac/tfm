 # rm(list = ls())

 source("ASWMain.R")

 
# Experimento 3 - Tratamiento con imputación de valores KNN, Drift, RSD y escalado de Pareto ####

object <- dataset

experimento3 <- object %>% 
  ASWImpute(method="knn") %>%
  ASWTransform(method="log", export=F) %>%
  ASWDriftCorrection(method="loess", plot=F, folderName="experimento3a/img/drifts/", export=F) %>%
  ASWTransform(method="power", export=F) %>%
  ASWFilter(method="filterByRSD", threshold=30, verbose=F) %>%
  ASWFilter(method="removeQCs") %>%
   ASWFilter(method="aggregateSamples") %>%
  ASWTransform(method="log") %>%
  ASWScale(method="pareto")

# saveRDS(experimento3, file = "experimento3a/rds/experimento3.rds")
# 
# experimento3 %>%
#   ASWPca(plot=T, scale=T, center=T)
# 
# 
# experimento3 %>%
#   ASWPca(plot=T)


experimento3 %>%
  ASWGenerateResults(experimentName="experimento3")
                     