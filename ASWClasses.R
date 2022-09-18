#' Classes definitions
#'
#' @description Classes like ASWExperiment, ASWStep, ASWPipeline and several related functions
#'
#' @author David Aguilera Castro
#'

# Class: ASWExperimentOptions ####
setClass("ASWExperimentOptions",
         slots = list(ncores = "numeric",
                      verbose = "logical",
                      round = "numeric"),
         prototype = list(ncores = bigstatsr::nb_cores(),
                          verbose=F,
                          round=4),)

# Class: ASWExperiment ####
setClass("ASWExperiment",
         slots = list(id = "character",
                      sample = "character",
                      codigo = "character",
                      order = "numeric",
                      style = "character",
                      features = "tbl_df",
                      metadata = "list",
                      options = "ASWExperimentOptions"
         ),
         prototype = list(metadata = as.list(NA)),
         # contains = "ASWPreTreatment"
)

# Method: loadExperiment ####
setGeneric("loadASWExperiment", function(object, filePath) {
  standardGeneric("loadASWExperiment")
})

loadASWExperiment <- function(filePath) {
            
            dataset <- read_excel(filePath)
            
            object <- new("ASWExperiment", 
                          id = dataset$ID, 
                          sample = dataset$Sample,
                          codigo = dataset$Codigo,
                          order = dataset$Order,
                          style = dataset$Style,
                          features = dplyr::select(dataset,starts_with("Feat")),
                          options = new("ASWExperimentOptions",verbose = T))
            
            rm('loadASWExperiment', pos = ".GlobalEnv")
            
            return (object)
          }

print.ASWExperiment <- function(object, ...){
  cat(primeColor("\nASW Experiment data object -----"))
  cat("\n· Número de muestras: ", length(object@style))
  cat("\n· Número de características: ", length(object@features))
  cat("\n· Características:\n", names(object@features))
}

# Method: getFeaturesSubset ####
setGeneric("getFeaturesSubset", function(object) {
  standardGeneric("getFeaturesSubset")
})

setMethod("getFeaturesSubset", "ASWExperiment",
          function(object){
            nFeatures <- length(object@features)
            
            subset <- c(
              seq(from=1, to=round(nFeatures*0.04)), 
              seq(from=(round(nFeatures*0.04)+4), to=(round(nFeatures*0.4)), by=8), 
              seq(from=(round(nFeatures*0.4)), to=nFeatures, by=20))
            
            return (subset)
          }
)

# Method: loadExperiment ####
setGeneric("ASWplot", function(p,filePath) {
  standardGeneric("ASWplot")
})

setMethod("ASWplot", "ANY",
          function(p, filePath){
            
            folderName <- dirname(filePath)
            
            ifelse(!dir.exists(file.path(".", folderName)), dir.create(file.path(".", folderName)), FALSE)
            
            png(filePath, width = 6000, height=5600, res = 600, pointsize = 12)
            
            print(p)
            
            dev.off()
          })

# Class: ASWStepOptions ####
setClass("ASWStepOptions",
         slots = list(
           verbose = "logical",
           export = "logical",
           plot = "logical"
         ),
         prototype = list(
           verbose = F,
           export = F,
           plot = F
         ),
)

# Class: ASWStep ####
setClass("ASWStep",
         slots = list(
           index = "numeric",
           method = "character",
           options = "ASWStepOptions"
         ),
         prototype = list(
           options = new("ASWStepOptions")
         ),
)

# Class: ASWStepImpute ####
setClass("ASWStepImpute",
         slots = list(
           threshold = "numeric"
         ),
         prototype = list(
           threshold = 70,
           index = 1,
           method="knn"
         ),
         contains = "ASWStep",
         validity = function(object) {
           if (!object@method %in% ASW@methodsAllowed$impute) {
             return("Method not allowed")
           }
         }
)

# Class: ASWStepTransform ####
setClass("ASWStepTransform",
         slots = list(),
         prototype = list(
           index = 2,
           method="log"
         ),
         contains = "ASWStep",
         validity = function(object) {
           if (!object@method %in% ASW@methodsAllowed$transform) {
             return("Method not allowed")
           }
         }
)

# Class: ASWStepDriftCorrection ####
setClass("ASWStepDriftCorrection",
         slots = list(),
         prototype = list(
           index = 3,
           method="loess"
         ),
         contains = "ASWStep",
         validity = function(object) {
           if (!object@method %in% ASW@methodsAllowed$drift) {
             return("Method not allowed")
           }
         }
         
)

# Class: ASWStepFilter ####
setClass("ASWStepFilter",
         slots = list(
           threshold = "numeric"
         ),
         prototype = list(
           threshold = 30,
           index = 4,
           method="filterByRSD"
         ),
         contains = "ASWStep",
         validity = function(object) {
           if (!object@method %in% ASW@methodsAllowed$filter) {
             return("Method not allowed")
           }
         }
         
)

# Class: ASWStepScale ####
setClass("ASWStepScale",
         slots = list(),
         prototype = list(
           index = 5,
           method="pareto"
         ),
         contains = "ASWStep",         
         validity = function(object) {
           if (!object@method %in% ASW@methodsAllowed$scale) {
             return("Method not allowed")
           }
         }
         
)

setGeneric("do", function(object, step) {
  standardGeneric("do")
})

setMethod("do", c("ASWExperiment", "ASWStepImpute"), 
          function(object, step){
            
            object <- ASWImpute(object,
                                method=step@method,
                                threshold=step@threshold,
                                verbose=step@options@verbose)

            if(step@options@export) {
              exportCSV(object,step)
            }

  return (object)
})

setMethod("do", c("ASWExperiment","ASWStepTransform"),
          function(object, step){
            object <- ASWTransform(object, 
                                   method=step@method)
  
            if(step@options@export) {
              exportCSV(object,step)
            }
  
          return (object)
})

setMethod("do", c("ASWExperiment", "ASWStepDriftCorrection"), 
          function(object,step){
            object <- ASWDriftCorrection(object, 
                                         method=step@method,
                                         verbose=step@options@verbose,
                                         plot=step@options@plot)
  
              if(step@options@export) {
                exportCSV(object,step)
              }
  
  return (object)
})

setMethod("do", c("ASWExperiment", "ASWStepFilter"), 
          function(object, step){
            object <- ASWFilter(object,
                                method=step@method,
                                threshold=step@threshold,
                                verbose=step@options@verbose)
  
            if(step@options@export) {
              exportCSV(object,step)
            }
  
            return (object)
          })

setMethod("do", c("ASWExperiment", "ASWStepScale"), 
          function(object, step){
            
          object <- ASWScale(object, 
                             method=step@method)
  
          if(step@options@export) {
            exportCSV(object,step)
          }
          
          return (object)
  
    })

# Method: exportCSV ####
setGeneric("exportCSV", function(object, step) {
  standardGeneric("exportCSV")
})

setMethod("exportCSV", c("ASWExperiment", "ASWStep"),
          function(object, step){
            
            # Folder
            folderName <- ifelse(missing(step), ASW@defaultOutputFolder , ASW@outputFolders[[step@index]])

            ifelse(!dir.exists(file.path(".", folderName)), dir.create(file.path(".", folderName)), FALSE)
            
            # File
            fileName <- ifelse(missing(step), ASW@defaultOutputFileName , step@method)
            filePath = paste(folderName, fileName, ".csv", sep="")
            
            matrix <- cbind(object@id, object@style, object@features)
            write.csv(matrix, filePath)
          }) 

# Class: ASWPipeline ####
setGeneric("ASWPipeline", function(object, steps) {
  standardGeneric("ASWPipeline")
})

setMethod("ASWPipeline", c("ASWExperiment", "list"),
          function(object, steps){
            
            for(i in 1:length(steps)){
              #Se ha pasado una lista de pasos en formato cadena, se intenta crear los pasos
              if(class(steps[[i]]) == "character") {
                stepClass = ASW@stepClasses[[names(steps)[i]]]
                stepMethod <- steps[[i]]
                
                step <- new(stepClass, method=stepMethod)
                
                object <- do(object,step)  
              }
              else {
                object <- do(object, steps[[i]])  
              }
            }
            
            return (object)
          })


# Method: ASWFeaturesIdentification ####
setGeneric("ASWFeaturesIdentification", function(object) {
  standardGeneric("ASWFeaturesIdentification")
})

setMethod("ASWFeaturesIdentification", "ASWExperiment",
          function(object){ 
            featuresID <-  read_excel(ASW@featuresIDFile)
            
            colnames(object@features) <- featuresID$id[featuresID$features %in% colnames(object@features)]
            
            return(object)
          
          })

# Method: ASWGenerateResults ####
setGeneric("ASWGenerateResults", function(object,experimentName, featuresSelectionMethods=c('none', 'splsda','rfe','fisher','vipR','vipS'), classifiers=c('plsda', 'lda', 'knn', 'rf'), nRepeats=1, fitTimes=10) {
  standardGeneric("ASWGenerateResults")
})

setMethod("ASWGenerateResults", "ASWExperiment",
          function(object, experimentName, featuresSelectionMethods, classifiers, nRepeats, fitTimes){
            #' @param experimentName "experimento2" 
            #' @param featuresSelectionMethods c('none', splsda','rfe','fisher','vipR','vipS')
            #' @param classifiers c('plsda', 'knn', 'rf')
            #' @param nRepeats 1
            #' @param fitTimes 10
            
            # Create folders ####
            experimentFolder <- paste(experimentName, "/", sep="")
            
            imgFolder   <- paste(experimentFolder, "img/", sep="")
            pcaFolder   <- paste(experimentFolder, "img/pca/", sep="")
            rdsFolder   <- paste(experimentFolder, "rds/", sep="")
            texFolder   <- paste(experimentFolder, "tex/", sep="")
            
            ifelse(!dir.exists(file.path(".", experimentFolder)), dir.create(file.path(".", experimentFolder)), FALSE)
            ifelse(!dir.exists(file.path(".", imgFolder)), dir.create(file.path(".", imgFolder)), FALSE)
            ifelse(!dir.exists(file.path(".", pcaFolder)), dir.create(file.path(".", pcaFolder)), FALSE)
            ifelse(!dir.exists(file.path(".", rdsFolder)), dir.create(file.path(".", rdsFolder)), FALSE)
            ifelse(!dir.exists(file.path(".", texFolder)), dir.create(file.path(".", texFolder)), FALSE)
            
            # PCA analysis ####
            object %>%
              ASWPca(plot=T, folderName=pcaFolder)
            
            # Features Selection ####
            
            fsFilePath <- paste(rdsFolder, experimentName, "fs.rds",sep="")
            
            if(file.exists(fsFilePath)) {
              fsObjects <- readRDS(fsFilePath)
              
              cat(primeColor("\n\n[FEATURES SELECTION] "))
              cat(red("done"))
            }
            else {
              fsObjects <- list()
              
              fsObjects$none <- object
              fsObjects$splsda <- ASWsPLSDAFeaturesSelection(object,verbose=T)
              fsObjects$rfe <- ASWRfeFeatureSelection(object,verbose=T)
              fsObjects$fisher <- ASWFFisherFeatureSelection(object,verbose=T)
              fsObjects$vipR <- ASWVIPbasedFeatureSelectionRelaxed(object,verbose=T)
              fsObjects$vipS <- ASWVIPbasedFeatureSelectionStrict(object,verbose=T)
              
              saveRDS(fsObjects, file = paste(rdsFolder, experimentName, "fs.rds",sep=""))
            }
            

            # Save FS results ####
            
            features <- data.frame(matrix(ncol = 2, nrow = 0))
            
            for(i in 1:length(fsObjects)) {
              fsMethod <- names(fsObjects)[i]
              featNumber <- length(fsObjects[[i]]@features)
              
              features <- rbind(features, c(fsMethod, featNumber))
            }
            
            features <- setNames(features, c("fsMethod", "features")) 
            
            print(xtable(features, type = "latex", include.rownames=F), file = paste(texFolder, experimentName, "fs.tex",sep=""))
            
            # Evaluation models ####
            resultadosBER <- data.frame(classifier=classifiers)
            resultadosAccuracy <- data.frame(classifier=classifiers)
            
            for(j in 1:length(featuresSelectionMethods)) {
              fsMethod <- featuresSelectionMethods[[j]]
              
              resultados = list()  
              
              cat(primeColor("\n\n[MODELS EVALUATION] "))
              cat(magenta(fsMethod))
              
              for(i in 1:length(classifiers)) {
                classifierName <- classifiers[[i]]
                
                cat(primeColor("\n\n[CLASSIFIER] "))
                cat(magenta(classifierName))
                
                resultado <- fsObjects[[fsMethod]] %>%
                  ASWTestRepeat(nrepeats=nRepeats, fitMethod=classifierName, fitTimes=fitTimes)
                
                resultados[[i]] <- resultado %>% purrr::list_modify("BERs" = NULL) %>% purrr::list_modify("Accuracies" = NULL)
              }
              
              resultadosDF = do.call(rbind.data.frame, resultados)
              resultadosBER[fsMethod] = resultadosDF$meanBER
              resultadosAccuracy[fsMethod] = resultadosDF$meanAccuracy
              
              saveRDS(resultadosDF, file = paste(rdsFolder,experimentName,"_",fsMethod,"_results.rds",sep=""))
              print(xtable(resultadosDF, type = "latex", include.rownames=FALSE), file = paste(texFolder,experimentName,"_",fsMethod,"_results.tex",sep=""))
            }
            
            saveRDS(resultadosBER, file = paste(rdsFolder,experimentName,"_results_BER.rds",sep=""))
            print(xtable(resultadosBER, type = "latex", include.rownames=FALSE), file = paste(texFolder,experimentName,"_results_BER.tex",sep=""))

            saveRDS(resultadosAccuracy, file = paste(rdsFolder,experimentName,"_results_Accuracy.rds",sep=""))
            print(xtable(resultadosAccuracy, type = "latex", include.rownames=FALSE), file = paste(texFolder,experimentName,"_results_Accuracy.tex",sep=""))
            
            }) 
