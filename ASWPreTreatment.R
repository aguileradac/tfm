#' Methods for pretreatment
#'
#' @description Methods for pretreatment steps.
#'
#' @author David Aguilera Castro


# Class: ASWPreTreatment ####
setClass("ASWPreTreatment",
         slots = list()
)

# Method: ASWExportCSV ####
setGeneric("ASWExportCSV", function(object, folderName, fileName) {
  standardGeneric("ASWExportCSV")
})

setMethod("ASWExportCSV", "ASWExperiment",
          function(object, folderName, fileName){
            
            ifelse(!dir.exists(file.path(".", folderName)), dir.create(file.path(".", folderName)), FALSE)
            
            filePath = paste(folderName, fileName, ".csv", sep="")
            
            matrix <- cbind(object@id, object@style, object@features)
            write.csv(matrix, filePath)
          }) 

# Method: ASWImpute ####
setGeneric("ASWImpute", function(object, ...) {
  standardGeneric("ASWImpute")
})

setMethod("ASWImpute", "ASWExperiment",
          function(object, method="KNN", threshold=70, verbose=F, export=F){
            #' @param method Impute method. 
            #' @param threshold Treshold for percentage of missing values allowed. Default: 70. 
            #' @param verbose Logical that indicates if more information by console should to be showed. Default: F.
            #' @param export Logical that indicates if the data matrix is to be exported to a csv file. Default: F.

            # Negative and Zero values ==> NA #### 
            object@features[object@features <= 0] <- NA
            
            # Remove features that exceed threshold of NA #### 
            features_to_supress <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("feat", "percent_na"))
            
            features_na = apply(X = is.na(object@features)*100, MARGIN = 2, FUN = mean)
            
            features_to_supress <- features_na[features_na > threshold]
            
            object@features <- dplyr::select(object@features, -names(features_to_supress))
            
            if (verbose) {
              cat(primeColor("\n\n[MISSING VALUE IMPUTATION] "))
              cat(magenta(method))
              
              cat(primeColor("\n[1] Features with a high percentage of missing values\n"))
              cat("\nThe following features shall be deleted for exceeding the threshold ")
              cat(red("[Threshold=",threshold,"%]\n", sep = ""))
              print(features_to_supress)
            }
            
            
            # Remove features that exceed threshold of NA by QC  #### 
            featuresByQC <- object@features[object@sample == "QC",]
            
            features_byqc_to_supress <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("feat", "percent_na"))
            
            features_byqc_na = apply(X = is.na(featuresByQC)*100, MARGIN = 2, FUN = mean)
            
            features_byqc_to_supress <- features_byqc_na[features_byqc_na > threshold]
            
            object@features <- dplyr::select(object@features, -names(features_byqc_to_supress))
            
            if (verbose) {
              cat(primeColor("\n[2] Features with a high percentage of missing values in QC\n"))
              cat("\nThe following features shall be deleted for exceeding the threshold ")
              cat(red("[Threshold=",threshold,"%]\n", sep = ""))
              print(features_byqc_to_supress)
              
              cat(primeColor("\n[3] Imputation algorithm\n"))
            }

            # Imputation methods ####
            if (method == ASW@methodsAllowed$impute$KNN) {
              imputedData <- impute::impute.knn(as.matrix(object@features))$data
            } 
            else if (method == ASW@methodsAllowed$impute$MF) {
              registerDoParallel(cores=detectCores()-1)
              imputedData <- missForest::missForest(as.matrix(object@features),verbose=verbose, parallelize='forest')$ximp
            } 
            else if (method == ASW@methodsAllowed$impute$LLS) {
              lls_estimation <- pcaMethods::llsImpute(Matrix = object@features, k = 10, center = TRUE, completeObs = TRUE, correlation = "pearson", allVariables = TRUE, maxSteps = 100)
              
              imputedData <- pcaMethods::completeObs(lls_estimation)
            }
            else if (method == ASW@methodsAllowed$impute$MEDIAN){
              imputedData <- apply(as.matrix(object@features), 2, function(x) {
                if(is.numeric(x)) ifelse(is.na(x), median(x, na.rm = TRUE),x) else x})
            }
            else if (method == ASW@methodsAllowed$impute$MEAN){
              imputedData <- apply(as.matrix(object@features), 2, function(x) {
                if(is.numeric(x)) ifelse(is.na(x), mean(x, na.rm = TRUE),x) else x})
            }
            else if (method == "remove") {
              featuresWithNA <- apply(X = is.na(object@features), MARGIN = 2, FUN = sum) 
              
              imputedData <- object@features[,(featuresWithNA == 0)]
            } 
            else if(method == "none") {
              imputedData <- object@features
            }
            
            # Return ####
            object@metadata$imputeMethod <- method
            object@features <- imputedData %>%
                                round() %>%
                                as_tibble()
            
            # CSV Export ####
            if(export) {
              ASWExportCSV(object, "output/impute/", method)
            }
            
            return (object)
          } 
)

# Method: ASWTransform ####
setGeneric("ASWTransform", function(object, ...) {
  standardGeneric("ASWTransform")
})

setMethod("ASWTransform", "ASWExperiment",
          function(object, method="log", export=F){
            #' @param method Transform method. 
            #' @param export Logical that indicates if the data matrix is to be exported to a csv file. Default: F.
            
            # Imputation methods ####
            if (method == ASW@methodsAllowed$transform$LOG) {
              transformFunction <- function(x) log(x)
            } 
            else if (method == ASW@methodsAllowed$transform$LOG2){
              transformFunction <- function(x) log2(x)
            }
            else if (method == ASW@methodsAllowed$transform$POWER){
              transformFunction <- function(x) round(exp(x),0)
            }
            else {
              transformFunction <- function(x) x
            }
            
            processedData <- apply(object@features, 2, transformFunction)      
            
            # Return ####
            object@features <- processedData %>% as_tibble()
            
            object@metadata$transformMethod <- method

            # CSV Export ####
            if(export) {
              ASWExportCSV(object, "output/transform/", method)
            }
            
            return (object)
          } 
)

# Method: ASWDriftCorrection ####
setGeneric("ASWDriftCorrection", function(object, ...) {
  standardGeneric("ASWDriftCorrection")
})

setMethod("ASWDriftCorrection", "ASWExperiment",
          function(object, method="loess", verbose=F, export=F, plot=F, folderName="images/drift/"){
            #' @param method Drift correction method. 
            #' @param verbose Logical that indicates if more information by console should to be showed. Default: F.
            #' @param export Logical that indicates if the data matrix is to be exported to a csv file. Default: F.
            #' @param plot Logical that indicates if graphs are to be produced. Default: F.
            #' @param folderName Path to save the graphs.
            
            if (verbose) {
              cat(primeColor("\n\n[DRIFT CORRECTION] "))
              cat(magenta(method))
            }
            
            ##
            
            # Build dataframe ####
            dsDrift = cbind("Sample" = object@sample, "Order" = object@order, object@features)
            rownames(dsDrift) <- object@id
            
            dsDriftCorrected = cbind("Sample" = object@sample, "Order" = object@order, object@features)
            rownames(dsDriftCorrected) <- object@id
            
            # Se obtiene los nombres de las filas de los QCs
            QCs <- object@id %>%
              as_tibble() %>%
              dplyr::filter(str_detect(value, "QC")) %>%
              pull(value)
            
            # Se obtienen los nombres de las características
            featNames <- names(object@features)

            meanQCsByFeature <- list()
            predictedY  <- list()
            predictedYC <- list()
            
            # Cálculo de Drift Correction ####
            for (i in 1:length(featNames)) {
              feat <- featNames[i]
              
              # Se calcula la media de los QCs para cada compuesto
              meanQCsByFeature[[i]] <- mean(dsDrift[QCs,feat])
              
              if(method == "loess") {
                # Se genera el modelo para cada característica
                model_loess <- stats::loess(formula = dsDrift[QCs,feat] ~ dsDrift[QCs,'Order'],data = dsDrift[QCs,], span = 0.5)
                
                # Se obtienen las predicciones
                predictedY[[i]] <- stats::predict(model_loess, dsDrift[,'Order'])
              } 
              else if(method == "SplineSmooth") {
                # Se genera el modelo para cada característica
                model_spline <- smooth.spline(dsDrift[QCs,'Order'], dsDrift[QCs,feat], spar=0.5)
                
                # Se obtienen las predicciones
                predictedY[[i]] <- stats::predict(model_spline, dsDrift[,'Order'])$y
              }
              
              # xcorrected(i) = xoriginal(i) + mean(Xqc) - fdrift(i)
              dsDriftCorrected[,feat] <- round(dsDrift[,feat] + meanQCsByFeature[[i]] - predictedY[[i]],object@options@round)
              
              if(method == "SplineSmooth") {
                model_splineC <- smooth.spline(dsDriftCorrected[QCs,'Order'], dsDriftCorrected[QCs,feat], spar=0.5)
                predictedYC[[i]] <- stats::predict(model_splineC, dsDriftCorrected[,'Order'])$y
                
              }
            }
            
            
            # Generación de gráficas
            if (plot) {
              cl = parallel::makeCluster(object@options@ncores)
              doParallel::registerDoParallel(cl)
              
              foreach (i=1:length(featNames), .packages=c('ggplot2','gridExtra')) %dopar% {
                if(i <= length(featNames)) {
                  feat <- featNames[i]
                  
                  # Gráfica antes del ajuste
                  fileName <- paste(folderName,feat,"_before.png", sep="")
                  
                  if(method == "loess") {
                  
                    before <- ggplot(dsDrift, aes(Order, get(feat))) +
                      geom_point(aes(colour = Sample)) +
                      geom_hline(aes(lty="mean", yintercept=meanQCsByFeature[[i]]), color = "green") +
                      geom_smooth(data = dsDrift[QCs,], method = method, mapping = aes_string(x='Order', y=feat))+
                      ylab(feat)
                  
                  } 
                  else if(method == "SplineSmooth") {
                    
                    before <- ggplot(dsDrift, aes(Order, get(feat))) +
                        geom_point(aes(colour = Sample)) +
                        geom_hline(aes(lty="mean", yintercept=meanQCsByFeature[[i]]), color = "green") +
                        geom_line(aes(x=dsDrift[,'Order'], y=predictedY[[i]]), color="steelblue2", size=1.25) +
                        ylab(feat)
                  }

                  
                  png(fileName, width = 6000, height = 5600, res = 600, pointsize = 12)
                  print(before)
                  dev.off()
                  
                  # Gráfica tras ajuste
                  fileName <- paste(folderName, feat,"_after.png", sep="")
                  
                  if(method =="loess") {
                    
                    after <- ggplot(dsDriftCorrected, aes(Order, get(feat))) + 
                      geom_point(aes(colour = Sample)) +
                      # geom_hline(aes(lty="mean", yintercept=meanQCs), color = "green") +
                      geom_smooth(data = dsDriftCorrected[QCs,], method = method, mapping = aes_string(x='Order', y=feat)) +
                      ylab(feat)
                    
                  } 
                  else if(method == "SplineSmooth") {
                    
                    after <- ggplot(dsDriftCorrected, aes(Order, get(feat))) + 
                      geom_point(aes(colour = Sample)) +
                      geom_line(aes(x=dsDriftCorrected[,'Order'], y=predictedYC[[i]]), color="steelblue2", size=1.25) +
                      ylab(feat)
                    
                  }                  
                  
                  png(fileName, width = 6000, height = 5600, res = 600, pointsize = 12)
                  print(after)
                  dev.off()
                  
                  # Gráfica combinada
                  fileName <- paste(folderName, feat,".png", sep="")
                  
                  png(fileName, width = 6000, height = 5600, res = 600, pointsize = 12)
                  print(grid.arrange(before, after, ncol=2))
                  dev.off()
                }
              }
              
              parallel::stopCluster(cl)
            }            
            
            # Return ####
            object@features <- dsDriftCorrected %>%
                               dplyr::select(starts_with("Feat")) %>% 
                               as_tibble()
            
            object@metadata$driftMethod <- method
            
            # CSV Export ####
            if(export) {
              ASWExportCSV(object, "output/drift/", method)
            }
            
            return (object)
          })

# Method: ASWFilter ####
setGeneric("ASWFilter", function(object, ...) {
  standardGeneric("ASWFilter")
})

setGeneric("ASWFilterByRSD", function(object, ...) {
  standardGeneric("ASWFilterByRSD")
})

setMethod("ASWFilterByRSD", "ASWExperiment", 
          function(object, threshold, verbose) {
            #' @param threshold Treshold for percentage of RSD allowed.
            #' @param verbose Logical that indicates if more information by console should to be showed.

            features_to_supress <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("feat", "percent_rsd"))
            
            # Datos de QC
            dfQC <- object@features[object@sample == "QC",]
            
            # Cálculo del coeficiente de variación para cada características
            rsd <- lapply(dfQC, function(x) (sd(x)/mean(x))*100)
            
            # Selección de características a eliminar
            features_to_supress <- rsd[rsd > threshold]
            
            if(verbose){
              print( rsd[rsd > threshold])
            }
            
            # Return ####
            object@features <- dplyr::select(object@features, -names(features_to_supress))
            
            return(object)
})

setGeneric("ASWFilterByRemoveQCs", function(object) {
  standardGeneric("ASWFilterByRemoveQCs")
})

setMethod("ASWFilterByRemoveQCs", "ASWExperiment", 
          function(object) {
            
            # Build return ####
            object@id <- object@id[object@sample != "QC"]
            object@codigo <- object@codigo[object@sample != "QC"]
            object@order <- object@order[object@sample != "QC"]
            object@style <- object@style[object@sample != "QC"]
            
            object@features <- object@features[object@sample != "QC",]
            object@sample <- object@sample[object@sample != "QC"]
            
            return(object)
          })


setMethod("ASWFilter", "ASWExperiment",
          function(object, method="filterByRSD", threshold=30, verbose=F, export=F){
            #' @param method Filter method. 
            #' @param threshold Treshold for percentage of RSD allowed. Default: 30. 
            #' @param verbose Logical that indicates if more information by console should to be showed. Default: F.
            #' @param export Logical that indicates if the data matrix is to be exported to a csv file. Default: F.
            
            if (verbose) {
              cat(primeColor("\n\n[FILTERING] "))
              cat(magenta(method))
            }
            
            ##
              
            if (method == ASW@methodsAllowed$filter$filterByRSD) {
              
              object <- ASWFilterByRSD(object, threshold, verbose)
              
            } else if (method == ASW@methodsAllowed$filter$removeQCs){
              
              object <- ASWFilterByRemoveQCs(object)

            } else if (method == ASW@methodsAllowed$filter$aggregateSamples){
              
              # Features
              features <- data.frame(matrix(ncol = length(object@features), nrow = 0))
              colnames(features) <- names(object@features)
              
              # Other arrays
              samples <- array()
              styles <- array()
              
              # Codigos
              codigos <- unique(object@codigo)
              
              for (i in 1:length(codigos)) {
                features[i,] = apply(object@features[object@codigo == codigos[i],],2,mean)
                samples[i] = object@sample[object@codigo == codigos[i]][1]
                styles[i] = object@style[object@codigo == codigos[i]][1]
              }
              
              # Return ####
              
              object@id <- codigos
              object@sample <- samples
              object@codigo <- codigos
              object@style <- styles
              object@order <- seq(from=1,to=length(codigos))
              object@features <- features %>% as_tibble()
            }

            # CSV Export ####
            if(export) {
              ASWExportCSV(object, "output/filter/", method)
            }
            
            return (object)
          } 
)


# Method: ASWScale ####
setGeneric("ASWScale", function(object, ...) {
  standardGeneric("ASWScale")
})

setMethod("ASWScale", "ASWExperiment",
          function(object, method="pareto", export=F){
            #' @param method Scale method. 
            #' @param export Logical that indicates if the data matrix is to be exported to a csv file. Default: F.

            # Scale methods ####
            if (method == ASW@methodsAllowed$scale$PARETO) {
              scaleFunction <- function(x) (x-mean(x))/sqrt(sd(x))
            } 
            else if (method == ASW@methodsAllowed$scale$AUTO){
              scaleFunction <- function(x) (x-mean(x))/sd(x)
            }
            else if (method == ASW@methodsAllowed$scale$RANGE){
              scaleFunction <- function(x) (x-mean(x))/(max(x) - min(x))
            }
            else if (method == ASW@methodsAllowed$scale$VAST){
              scaleFunction <- function(x) ((x-mean(x))/sd(x))*(mean(x)/sd(x))
            }
            else if (method == ASW@methodsAllowed$scale$LEVEL){
              scaleFunction <- function(x) ((x-mean(x))/sd(x))*(mean(x)/mean(x))
            }
            
            processedData <- round(apply(object@features, 2, scaleFunction),object@options@round)       
            
            # Return ####
            object@features <- processedData %>% as_tibble()
            
            object@metadata$scaleMethod <- method
            
            
            # CSV Export ####
            if(export) {
              ASWExportCSV(object, "output/scale/", method)
            }
            
            return (object)
          } 
)

# Method: ASWPca ####
setGeneric("ASWPca", function(object, ...) {
  standardGeneric("ASWPca")
})

setMethod("ASWPca", "ASWExperiment",
          function(object, ncomp=5, center=F, scale=F, plot=F, folderName = 'images/pca/'){
            #' @param ncomp Number of components for PCA. Default: 5. 
            #' @param center Logical that indicates if data will be centered. Default: F.
            #' @param scale Logical that indicates if data will be scaled. Default: F.
            #' @param plot Logical that indicates if graphs are to be produced. Default: F.
            #' @param folderName Path to save the graphs.
            
            X <- as.matrix(object@features)
            Y <- as.factor(object@style)
            
            pca <- mixOmics::pca(X, ncomp = ncomp, center = center, scale = scale)
            
            cat(primeColor("\n\n[PRETREATMENT EVALUATION] "))
            cat(magenta("PCA"))
            
            imputeMethod <- ifelse(is.null(object@metadata$imputeMethod), "none", object@metadata$imputeMethod)
            transformMethod <- ifelse(is.null(object@metadata$transformMethod), "none", object@metadata$transformMethod)
            scaleMethod <- ifelse(is.null(object@metadata$scaleMethod), "none", object@metadata$scaleMethod)
            driftMethod <- ifelse(is.null(object@metadata$driftMethod), "none", object@metadata$driftMethod)
            
            cat(paste("\nImputation method:",imputeMethod))
            cat(paste("\nDrift correction method:",driftMethod))
            cat(paste("\nTransformation method:",transformMethod))
            cat(paste("\nScale method:",scaleMethod))
            # cat(paste("\nDrift correction method:",dc))
            # cat(paste("\nRSD cutoff:",rf))
            cat("\n--------------------------------\n")
            
            print(pca$prop_expl_var$X)
            
            # Plot PCA ####
            if(plot) {
              col.pch <- list(col = c("dodgerblue", "firebrick", "forestgreen", "gold"), pch = c())
              IDs = object@id
              
              ifelse(!dir.exists(file.path(".", folderName)), dir.create(file.path(".", folderName)), FALSE)
              
              ##
              
              # Se generan las gráficas
              cl = parallel::makeCluster(object@options@ncores)
              doParallel::registerDoParallel(cl)
              
              foreach (i=2:ncomp, .packages='mixOmics') %dopar% {
              # for (i in 2:ncomp) {
                comp = c(1,i)
                
                # filenameScores = paste(folderName, "scores1_",i,".tiff", sep="")
                # filenameLoadings = paste(folderName, "loadings1_",i,".tiff", sep="")
                filenameScores = paste(folderName, "scores1_",i,".png", sep="")
                filenameLoadings = paste(folderName, "loadings1_",i,".png", sep="")
                
                # Scores
                scoresPlot <- mixOmics::plotIndiv(pca, group = Y, ind.names = IDs, legend = T, cex = 1, style = "graphics",
                                    abline = TRUE, title = "", size.xlabel = 1.7, comp = comp, size.axis = 1.2,size.legend = 1,
                                    size.title = 2, star = F, col = col.pch$col[1:length(levels(Y))])

                png(filenameScores, width = 6000, height=5600, res = 600, pointsize = 12)
                #tiff(filenameScores, width = 6000, height = 5600, res = 600, compression = "lzw", pointsize = 12)
                mixOmics::plotIndiv(pca, group = Y, ind.names = IDs, legend = T, cex = 1, style = "graphics",
                                    abline = TRUE, title = "", size.xlabel = 1.7, comp = comp, size.axis = 1.2,size.legend = 1,
                                    size.title = 2, star = F, col = col.pch$col[1:length(levels(Y))])
                
                # mtext(paste("Imputation method:",imputeMethod,
                #             "Drift correction method:",driftMethod,
                #             "Transformation method:",transformMethod,
                #             "Scale method",scaleMethod), side=3)
                dev.off()
                
                # Loadings
                mixOmics::plotVar(pca, cex = 1, style = "graphics", abline = TRUE, title = "", comp = comp, var.names = F, legend = T)
                
                png(filenameLoadings, width = 6000, height=5600, res = 600, pointsize = 12)
                #tiff(filenameLoadings, width = 5400, height = 5300, res = 600, compression = "lzw", pointsize = 13)
                mixOmics::plotVar(pca, cex = 1, style = "graphics", abline = TRUE, title = "", comp = comp, var.names = F, legend = T)
                dev.off()
                
              }
              
              parallel::stopCluster(cl)
              
            }
          } 
)

