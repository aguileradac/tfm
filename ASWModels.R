#' Models training and testing
#'
#' @description Methods for fit and test classification models
#'
#' @author David Aguilera Castro
#'

# Method: ASWFitModel ####
setGeneric("ASWFitModel", function(object, ...) {
  standardGeneric("ASWFitModel")
})

setMethod("ASWFitModel", "ASWExperiment",
          function(object, method="plsda", partition=0.7, times=1, plot=F, folderName="images/plsda/"){
            #' @param method Classification model. Options are: "plsda", "knn", "rf" or "lda". 
            #' @param partition Size of data partition for training model. Default: 0.7. 
            #' @param times Number of partitions. Default: 1.
            #' @param plot Logical that indicates if graphs are to be produced. Default: F.
            #' @param foldeName Path to save the graphs. Default: "images/plsda/".
            
            X <- object@features
            Y <- as.factor(object@style)
            
            trainIndex <- createDataPartition(Y, p = partition, list = T, times = times)
            
            Xtrain <- list()
            Xtest  <- list()
            Ytrain <- list()
            Ytest  <- list()
            
            cat(primeColor(" Training: "))
            
            for (i in 1:length(trainIndex)) {
              currentTrainIndex <- unlist(trainIndex[i])
              
              Xtrain[[i]] <- X[currentTrainIndex,]
              Xtest[[i]]  <- X[-currentTrainIndex,]
              
              Ytrain[[i]] <- Y[currentTrainIndex]
              Ytest[[i]]  <- Y[-currentTrainIndex]
            }
            
            # Train control ####
            ctrl <- trainControl(method="LOOCV", allowParallel = F, classProbs = T, summaryFunction = multiClassSummary, verboseIter = F) #,classProbs=TRUE,summaryFunction = twoClassSummary)
            
            # Hyperparameter tuning ####
            
            if(method == "knn") {
              
              tuneGrid <- expand.grid(k = seq(1,10))
              
            } else if (method == "lda") {
              
              tuneGrid <- NULL
              
            } else if (method == "pls") {
              
              tuneGrid <- expand.grid(ncomp = seq(1,5))
            
            } else if (method == "plsda") {
              
              ncomp = 5
              
            } else if (method == "rf") {
              
              mtry <- sqrt(ncol(Xtrain[[i]]))
              
              tuneGrid <- expand.grid(mtry=1:mtry)

            } else if (method == "svmLinear") {
              
              tuneGrid <- expand.grid(C = c(0.25, 0.5, 1, 2, 4))
              
            } else if (method == "svmLinear2") {
              
              tuneGrid <- expand.grid(cost = seq(0, 4, by=0.5))
              
            } else if (method == "svmPoly") {
              
              tuneGrid <- expand.grid(scale = c(0.1, 0.2),C = c(0.75, 1, 1.25),degree = 0:4)
              
            } else {
              
              tuneGrid <- NULL  
            }
            
            # Model training ####
            modelFit <- list()
            
            for(i in 1:length(trainIndex)) {
              cat("*")

              
              if(method == 'plsda') {
                  result.plsda <- mixOmics::plsda(Xtrain[[i]], Ytrain[[i]], ncomp=ncomp) # run the method
                
                  perf.plsda <- perf(result.plsda, validation = "Mfold", nrepeat = 10, folds = 6,
                                    progressBar = F, auc = TRUE) # include AUC values
                
                  optimalNComp <- perf.plsda$choice.ncomp[2,3]
                  
                
                  modelFit[[i]] <- mixOmics::plsda(Xtrain[[i]], Ytrain[[i]], ncomp=optimalNComp) # run the method
                  
                  if(plot) {
                    p1 <- plotIndiv(modelFit[[i]], ind.names = FALSE, legend=TRUE,
                              ellipse = F, star = TRUE, title = 'PLS-DA Model')
                    
                    ASWplot(p1, 
                            paste(folderName, "train_",i,".png", sep=""))
                    
                    #
                    
                    fileName = paste(folderName, "loadings_",i,".png", sep="")

                    png(fileName, width = 6000, height=6000, res = 600, pointsize = 12)
                    plotVar(modelFit[[i]], legend=TRUE,title = 'PLS-DA Model', X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
                    dev.off()
                    
                    #
                    
                    background <- background.predict(modelFit[[i]], comp.predicted=2,dist = "max.dist") 
                    
                    p3 <-  plotIndiv(modelFit[[i]], comp = 1:2, group = Ytrain[[i]], ind.names = FALSE, 
                                    title = "Maximum distance", legend = TRUE,  background = background)
                    
                    ASWplot(p3, 
                            paste(folderName, "background_",i,".png", sep=""))
                    
                    #
                    
                    legend=list(legend = levels(Ytrain[[i]]), # set of classes
                                col = unique(color.mixo(Ytrain[[i]])), # set of colours
                                title = "Wine Type", # legend title
                                cex = 0.7) # legend size
                    
                    fileName = paste(folderName, "heatmap_",i,".png", sep="")
                    
                    png(fileName, width = 6000, height=5600, res = 600, pointsize = 12)
                    cim(modelFit[[i]], row.sideColors = color.mixo(Ytrain[[i]]), legend = legend)
                    dev.off()
                    
                  }
                

              } else {
                modelFit[[i]] <- caret::train(x=Xtrain[[i]], y=Ytrain[[i]], 
                                            method = method, 
                                            trControl = ctrl, 
                                            tuneGrid = tuneGrid,
                                            metric = 'Accuracy')
              
              }
            }
            
            ## RESULTS
            
            models <- list()
            
            models$fit <- modelFit
            models$X <- X
            models$Y <- Y
            models$Xtrain <- Xtrain
            models$Ytrain <- Ytrain
            models$Xtest <- Xtest
            models$Ytest <- Ytest
            
            return(models)
            
          }
)

# Method: ASWTestModel ####
setGeneric("ASWTestModel", function(models, export=F) {
  standardGeneric("ASWTestModel")}
)

setMethod("ASWTestModel", c("list"),
          function(models, export){
            #' @param models List of trained models.  

            Ypredicted <- list()
            confusionMatrix <- list()
            classifier_metrics <- list()
            
            models.list <- list()
            results <- list()
            bestTunes <- list() 
            
            meanBER   <- 0
            
            for(i in 1:length(models$fit)) {
              model <- models$fit[[i]]
              
              if (is(model) == "mixo_plsda") {
                Ypredicted1 = predict(model, newdata = models$Xtest[[i]])
                Ypredicted[[i]] = Ypredicted1$class$mahalanobis.dist[,ncol(Ypredicted1$class$mahalanobis.dist)]
              } else {
                Ypredicted[[i]] <- predict(model, newdata = models$Xtest[[i]])  
              }
              

              confusionMatrix[[i]] = get.confusion_matrix(truth = models$Ytest[[i]], predicted = Ypredicted[[i]]) 
              classifier_metrics[[i]] <- ml_test(Ypredicted[[i]], models$Ytest[[i]], output.as.table = FALSE)
              
              result <- list()
              
              result$model <- model
              result$confusionMatrix <- confusionMatrix[[i]]
              result$classifier_metrics <- classifier_metrics[[i]]
              
              # BER
              BER <- get.BER(confusionMatrix[[i]])
              
              result$classifier_metrics$BER <- BER
              
              meanBER <- meanBER + BER
              
              result$Xtrain <- models$Xtrain[[i]]
              result$Ytrain <- models$Ytrain[[i]]
              result$Xtest <- models$Xtest[[i]]
              result$Ytest <- models$Ytest[[i]]
              
              models.list[[i]] <- result
              
              bestTunes[[i]] <- model$bestTune
            }
            
            results$X <- models$X
            results$Y <- models$Y
            
            results$models <- models.list
            
            listTunes <- lapply(bestTunes, c, recursive=TRUE)

            
            # Cálculo BER  
            # results$meanBER <- round(meanBER/length(models$fit),2)
            
            summariseMetrics <- results$models %>%
              flatten() %>%
              tibble::enframe() %>%
              filter(name == 'classifier_metrics') %>%
              unnest_wider(value) 
            
            results$meanBER <- summariseMetrics %>%
              summarise(mean = round(mean(BER),2)) %>%
              pull(mean)
            
            results$meanAccuracy <- summariseMetrics %>%
              summarise(mean = round(mean(accuracy),2)) %>%
              pull(mean)
            
            results$meanRecall <- summariseMetrics %>%
              dplyr::select(recall) %>%
              unnest_wider(recall) %>%
              dplyr::summarise(across(everything(), mean))
            
            results$meanSpecificity <- summariseMetrics %>%
              dplyr::select(specificity) %>%
              unnest_wider(specificity) %>%
              dplyr::summarise(across(everything(), mean))
            
            #####
            
            class(results) <- "ASWTestModelResults"
            
            #####
            
            return(results)
            
          }
)

print.ASWTestModelResults <- function(results, ...){
  cat("El parámetro óptimo medio es",results$optimumParameter,"\n")
  cat("mean(BER): ",results$meanBER,"\t")
  cat("mean(accuracy): ",results$meanAccuracy,"\n")
  cat("mean(recall): ","\n")
  print(results$meanRecall)
  cat("mean(specificity): ","\n")
  print(results$meanSpecificity)
}

# Method: ASWTestRepeat ####
setGeneric("ASWTestRepeat", function(object, nrepeats=1, fitMethod="pls", fitTimes=1) {
  standardGeneric("ASWTestRepeat")}
)

setMethod("ASWTestRepeat", "ASWExperiment",
          function(object, nrepeats, fitMethod, fitTimes){
            
          meanBERs <- array()
          meanAccuracies <- array()
          meanTime <- array()
          
          for(i in 1:nrepeats) {
            
              cat(primeColor("\nRepeat: "))
              cat(i,"/",nrepeats, sep="")
            
              tic = Sys.time()
              
              result<- object %>%
                ASWFitModel(method=fitMethod,times=fitTimes) %>%
                  ASWTestModel()
              
              toc = Sys.time()
              
              meanBERs[i] <- result$meanBER
              meanAccuracies[i] <- result$meanAccuracy
              meanTime[i] <- (toc-tic)
              
              cat(primeColor(" Time elapsed: "))
              cat(HMS(as.numeric(difftime(toc,tic,units='secs'))))
              
           }
              
            #####
            result <- list()  
          
            result$nrepeats   <- nrepeats
            result$fitMethod  <- fitMethod
            result$fitTimes   <- fitTimes
            result$meanTime      <- mean(meanTime)
            result$totalTime      <- sum(meanTime)
            
            result$BERs <- meanBERs
            result$meanBER      <- mean(meanBERs)
            

            result$Accuracies <- meanAccuracies
            result$meanAccuracy <- mean(meanAccuracies)
            
            class(result) <- "ASWTestRepeatResults"
            
            #####
            
            # beep()
            
            return(result)
            
          }
)

print.ASWTestRepeatResults <- function(result, ...){
  cat("\n\nTestRepeat Results -----")
  cat("\nNúmero de repeticiones: ",result$nrepeats)
  cat("\nFit Model: ",result$fitMethod,"\tFit times: ",result$fitTimes)
  cat("\nMean duration: ",result$meanTime,"\tTotal duration: ",result$totalTime)
  cat("\n------------------------")
  cat("\nBERs: ",head(result$BERs,10))
  cat("\nBER mean: ",result$meanBER)
  cat("\n\nAccuracy: ",head(result$Accuracies,10))
  cat("\nAccuracy mean: ",result$meanAccuracy)
  cat("\n------------------------\n")
}
