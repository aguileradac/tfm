#' Methods for feature selection
#'
#' @description Methods for feature selection step.
#'
#' @author David Aguilera Castro

# Method: ASWFeaturesSelection ####
setGeneric("ASWFeaturesSelection", function(object, method, verbose=F) {
  standardGeneric("ASWFeaturesSelection")
})

setMethod("ASWFeaturesSelection", "ASWExperiment",
          function(object, method, verbose){
            
            if(method == "ALL") {
              methods = c('rfe', 'vipS', 'vipR', 'splsda', 'fisher')
            } else {
              methods = c(method)
            }
            
            objects <- list()
            
            for(i in 1:length(methods)) {
              if (methods[i] == "splsda") {
                objects$splsda <- ASWsPLSDAFeaturesSelection(object, verbose)
              } else if (methods[i] == "rfe") {
                objects$rfe <- ASWRfeFeatureSelection(object, verbose)
              } else if (methods[i] == "fisher") {
                objects$fisher <- ASWFFisherFeatureSelection(object, verbose)
              } else if (methods[i] == "vipS") {
                objects$vipS <- ASWVIPbasedFeatureSelectionStrict(object, verbose)
              } else if (methods[i] == "vipR") {
                objects$vipR <- ASWVIPbasedFeatureSelectionRelaxed(object, verbose)
              }
              
              unregister_dopar()
            }
            
            return (objects)
            
          }) 

# Method: ASWsPLSDAFeaturesSelection ####
setGeneric("ASWsPLSDAFeaturesSelection", function(object, verbose=F, export=F, plot=F) {
  standardGeneric("ASWsPLSDAFeaturesSelection")
})

setMethod("ASWsPLSDAFeaturesSelection", "ASWExperiment",
          function(object, verbose, export, plot){
            
            if (verbose) {
              cat(primeColor("\n[FEATURES SELECTION] "))
              cat(magenta("sPLS-DA "))
            }
            
            ##  
            
            X <- as.matrix(object@features)
            rownames(X) <- object@codigo
            Y <- as.factor(object@style)

            components = 5
            validation = "Mfold"
            folds = 5
            nrepeat = 10

            labels = F
            ellipse = T
            legend_position = "bottom"
            
            # Pre-process ####
            
            list_keepX <- getFeaturesSubset(object)

            # Tuning ####
            
            tune_splsda <- mixOmics::tune.splsda(X, Y, ncomp = components, validation = validation, folds = folds,
                                                 progressBar = verbose, dist = 'max.dist', measure = "BER",
                                                 test.keepX = list_keepX, nrepeat = nrepeat, cpus = detectCores()-1) # cpus = 4
            
            error <- tune_splsda$error.rate
            ncomp <- tune_splsda$choice.ncomp$ncomp # optimal number of components based on t-tests
            select_keepX <- tune_splsda$choice.keepX[1:ncomp]  # optimal number of variables to select
            
            errors_splsda <- data.frame(tune_splsda$error.rate) %>% 
              rownames_to_column("feature") %>%
              pivot_longer(cols = -feature)
            
            errors_sd <- data.frame(tune_splsda$error.rate.sd) %>% 
              rownames_to_column("feature_sd") %>%
              pivot_longer(cols = -feature_sd)
            
            errors_splsda <- cbind(errors_splsda, sd = errors_sd$value) %>% 
              as.data.frame() %>% 
              as_tibble()
            
            bal_error_rate <- ggplot(data = errors_splsda, aes(x = feature, y = value, group = name)) +
              geom_line(aes(color = name)) +
              geom_point(aes(color = name)) +
              geom_errorbar(aes(ymin = value-sd, ymax = value+sd, color = name), width=.1) +
              theme_bw() +
              xlab("Number of features") +
              ylab("Error") +
              theme(legend.title = element_blank(),
                    legend.position = legend_position) +
              scale_colour_viridis_d(begin = 0, end = 0.8)
            
            # Final model ####
            
            if (ncomp == 1){
              ncompX <- 2
            }else{
              ncompX <- ncomp}
            
            ####
            
            res_splsda <- mixOmics::splsda(X, Y, ncomp = ncompX, keepX = select_keepX)
            
            SPLSDAi <- data.frame(res_splsda$variates$X, Groups = Y) %>% 
              rownames_to_column("ID")
            
            splsda_scores_plot <- ggplot(SPLSDAi, aes(x = comp1, y = comp2, color = Groups, shape = Groups, label = ID)) +
              {if(!labels)geom_point(size = 2, alpha = 0.9)} +
              xlab("Component 1") +
              ylab("Component 2") +
              {if(ellipse)stat_ellipse(type = "norm")} +
              {if(labels)geom_text(aes(label = ID), show.legend = TRUE)} +
              theme_bw() +
              theme(legend.title = element_blank(),
                    legend.position = legend_position) +
              scale_colour_viridis_d(begin = 0, end = 0.8)
            

            scores_splsda <- SPLSDAi %>% 
              dplyr::select(-Groups, -ID) %>%
              as_tibble()
            
            selected_variables <- mixOmics::selectVar(res_splsda, comp = 1)
            selected_variables <- round(selected_variables$value, 4)
            selected_variables <- data.frame(feature = rownames(selected_variables), 
                                             value = selected_variables$value.var) %>% 
              as_tibble()
            
            # Return ####
            
            object@features <- object@features[,colnames(object@features) %in% selected_variables$feature]
            
            object@metadata$featuresSelectionMethod <- "splsda"
            
            # Plot ####
            if(plot){
              print(bal_error_rate)
              print(splsda_scores_plot)
            }
            
            
            
            
            # CSV Export #### 
            if(export) {
              ASWExportCSV(object, "output/fs/", "splsda")  
            }   
            
            return (object)
          }
          
)

# Method: ASWRfeFeatureSelection ####
setGeneric("ASWRfeFeatureSelection", function(object, verbose=F, export = F, plot = F) {
  standardGeneric("ASWRfeFeatureSelection")
})

setMethod("ASWRfeFeatureSelection", "ASWExperiment",
          function(object, verbose, export, plot){

            if (verbose) {
              cat(primeColor("\n[FEATURES SELECTION] "))
              cat(magenta("Recursive Features Selection\n"))
            }
            
            ##  
            
            X <- as.matrix(object@features)
            rownames(X) <- object@codigo
            Y <- as.factor(object@style)
            
            # Control ####
            control <- rfeControl(functions=rfFuncs, method="cv", number=5, allowParallel=T)  
            
            # Process ####
            list_keepX <- getFeaturesSubset(object)
            
            results <- rfe(X, Y, sizes=list_keepX, rfeControl=control)
            
            # Return ####
            features_selected = predictors(results)  
            
            dfData <- object@features %>%
              dplyr::select(all_of(features_selected[order(features_selected)]))

            object@features <- dfData %>% as_tibble()

            object@metadata$featuresSelectionMethod <- "rfe"
            
            # CSV Export #### 
            if(export) {
              ASWExportCSV(object, "output/fs/", "rfe")  
            }   
            
            return (object)
          }
          
)

# Method: ASWVIPbasedFeatureSelection ####
setGeneric("ASWVIPbasedFeatureSelection", function(object, frequenceMinStability = 0.9, verbose=F, export=F, plot=F) {
  standardGeneric("ASWVIPbasedFeatureSelection")
})

setMethod("ASWVIPbasedFeatureSelection", "ASWExperiment",
          function(object, frequenceMinStability, verbose, export, plot){
            
            if (verbose) {
              cat(primeColor("\n[FEATURES SELECTION] "))
              cat(magenta("VIPbased Features Selection "))
              cat(red("[frequenceMinStability=",frequenceMinStability,"%]\n", sep = ""))
            }
            
            # Auxiliar functions ####
            
            iterSelect = function(X, Y, var.select, seq.iter, repetcv, compPLSDA){

              Ymodel = Y
              all.select.var = all.VIPfinal = all.Nvar.select = all.compModr = all.BER.r = c()
              
              for (i in 1:repetcv) {
                all.select.var = c(all.select.var, var.select[[i]]$select.Var1)
                all.VIPfinal = c(all.VIPfinal, var.select[[i]]$VIPfinal)
                all.Nvar.select = c(all.Nvar.select, var.select[[i]]$Nvar.select)
                all.compModr = c(all.compModr, var.select[[i]]$compModround)
                all.BER.r = c(all.BER.r, var.select[[i]]$BER.r)
              }
              
              compPLSDAini = compPLSDA
              variables=as.matrix(table(all.select.var))
              info.modelo=cbind(all.Nvar.select, all.VIPfinal, all.compModr, all.BER.r)
              colnames(info.modelo)=c("Compuestos R", "VIP", "Componentes R", "BER R")
              rownames(info.modelo) = rownames(info.modelo, do.NULL = FALSE, prefix = "Iter")
              selection.freq = matrix(-1,compPLSDA,length(seq.iter))
              
              #' Parallel VIP iteration
              
              cl = parallel::makeCluster(object@options@ncores)
              doParallel::registerDoParallel(cl)
              
              iterFreq = foreach(iter = 1:length(seq.iter), .packages = "mixOmics") %dopar% {
                Xiter = X
                Frequence_minStability = seq.iter[iter]
                minStability = Frequence_minStability * repetcv
                remainVar = c()
                choise=rownames(variables)[which(variables>minStability)]
                for (choisei in 1:length(choise)){
                  remainVar[choisei] = which(colnames(X)==choise[choisei])
                }
                
                if(length(remainVar) < 2){
                  compPLSDAini
                  selection.freq = as.matrix(rep(NA, compPLSDAini))
                }else{
                  Xiter = X[,remainVar]
                  
                  if (ncol(Xiter)<=compPLSDA) {
                    compPLSDA = ncol(Xiter) - 1
                  }
                  
                  my_plsda=plsda(Xiter,Ymodel, ncomp = compPLSDA, scale = TRUE)
                  perf_plsda=perf(my_plsda, dist = "mahalanobis", auc = TRUE, progressBar = FALSE, validation = "loo")
                  selection.freq = perf_plsda$error.rate$BER
                  selection.freq
                  rownames = rownames(perf_plsda$error.rate$BER)
                }
                
                return(list("selection.freq" = selection.freq, "rownames" = rownames))
                
              }
              
              parallel::stopCluster(cl)

              selection.freq = matrix(ncol = length(seq.iter), nrow = compPLSDA)
              for (iter in 1:length(seq.iter)) {
                selection.freq[1:nrow(iterFreq[[iter]]$selection.freq),iter] = iterFreq[[iter]]$selection.freq
              }
              
              colnames(selection.freq) = seq.iter
              rownames(selection.freq) = iterFreq[[1]]$rownames
              selection.freq
              select.error = selection.freq[selection.freq>=0]
              optimal.error = min(select.error[!is.na(select.error)])
              #########################################
              
              return(list("model.info" = info.modelo, "variables" = variables, "selection.freq" = selection.freq,
                          "optimal.error" = optimal.error))
              
            }
            
            choiseComp = function(choisecomp, minchoisecomp){
              if (!exists("minchoisecomp")) {
                minchoisecomp = 0.05
              }
              if (min(choisecomp) > minchoisecomp) {
                ncomp = which(choisecomp == min(choisecomp))[1]
              }else{
                ncomp = which(choisecomp <= minchoisecomp)[1] #Elijo componentes por el error balanceado de la mahalanobis.dist
              }
              if(length(ncomp)>1){ncomp=ncomp[1]}
              
              #Numero de componentes en cada iteracion
              return(ncomp)
            }
            
            # VIPbased FS ####
            
            X <- object@features
            Y <- object@style
            
            # Parameters
            repetcv = 30 # Number of repetitions in the outer loop cross-validation 
            compPLSDA = 5 # Number of components of the PLS-DA
            nfold = 6 # (8) n-fold outer loop cross-validation 
            BERc = 0.8 # Cutoff BER. See mean overall BER in double cross validation
            VIPinic = 0.8 # Starting iteration
            saltoVIP = 0.05 # Increasing of VIP in each iteration
            oportunidad = 1 # Number of chances. Select 1 or 2
            minimVar = 10 # Minimal number of variables selected
            seq.iter = seq(0.1,0.98,0.02) # 

            selectVIP = c()
            Ymodel = Y
            
            
            message("Obtención de var.select")
            
            # Obtención de var.select ####
            {
              var.select <- list()
              
              # cl = parallel::makeCluster(object@options@ncores)
              # doParallel::registerDoParallel(cl)
              # dcvalidation = foreach(cv = 1:repetcv, .packages = "mixOmics") %dopar% {
              # foreach(cv = 1:repetcv, .packages = "mixOmics") %dopar% {
              for(cv in 1:repetcv) {
                if(verbose) {
                  flush.console()
                  cat("[CV] ",cv,"/",repetcv," \r",sep="")
                  # print(paste("CV: ", cv,sep=""))
                }
                #####################################################################################
                #Particion del test set. El rest set son el resto
                #####################################################################################
                level=levels(as.factor(Ymodel)) 
                test=c()
                nfactors=c()
                fold = 1
                #Elijo training y test set
                for (niveles in 1:length(level)){
                  factors=round(length(which(Ymodel == level[niveles]))/nfold)
                  #Elijo todos los test y en el ultimo elijo el resto
                  if (fold==nfold){
                    totalfactor=length(which(Ymodel == level[niveles]))
                    test = c(test, sample(which(Ymodel == level[niveles]))[c(totalfactor - (totalfactor-c(factors*(nfold-1)+1))): totalfactor])
                  }else{
                    test = c(test, sample(which(Ymodel == level[niveles]))[c(factors*fold-factors+1):c(factors*fold)])
                  }
                }
                train = c(1:length(Ymodel))[-test] #El train set es el resto que no esta en el test set
                
                #Construyo el primer modelo y cojo los VIPs
                my_plsda_train = mixOmics::plsda(X[train,], Ymodel[train], ncomp = 4)
                my_vip=vip(my_plsda_train) 
                
                for(selC in 1:dim(my_vip)[1]){
                  selectVIP[selC]=my_vip[selC,which.max(my_vip[selC,])]
                } 
                as.matrix(selectVIP) #Primero ponemos todos los VIP mas altos en alguna LV en una misma matriz
                
                #Ahora comienza la iteracion con los diferentes VIP mientras mejore mi modelo 
                VIPinicial=VIPinic
                continuar=0
                for(threshold in c(seq(VIPinicial,2,saltoVIP))){
                  select.Var=which(selectVIP>threshold) #Seleccion de variables del modelo inicial con todas las variables
                  if (length(select.Var)<minimVar){break} #Parar si tengo menos de X compuestos elegidos (overfitting)
                  if (threshold==VIPinicial){
                    Xselect=t(t(X)[select.Var,])
                  } else{
                    Xselect=t(t(Xselect)[select.Var,])
                  }
                  compMod=c()
                  BER.select=c()
                  
                  if (threshold == VIPinicial){
                    BER = BERc
                  }
                  
                  Xsel=Xselect[train,] #Selecciono las muestras
                  Ysel=Y[train] #Selecciono las muestras
                  
                  my.plsda.select <- mixOmics::plsda(Xsel, Ysel, ncomp = compPLSDA)
                  perf.plsda.select <- perf(my.plsda.select, validation = "loo", progressBar = FALSE)
                  choisecomp=as.data.frame(perf.plsda.select$error.rate$BER)[,1]
                  # source("./mane/choiseComp.R")
                  ncomp = choiseComp(choisecomp = choisecomp, minchoisecomp = 0.03)
                  
                  my.plsda.perf.select=mixOmics::plsda(Xsel, Ysel, ncomp = ncomp)
                  
                  test.predict.select <- predict(my.plsda.perf.select, Xselect[test,], dist = "mahalanobis.dist") #Prediccion del modelo
                  prediction.select = test.predict.select$class$mahalanobis.dist[,ncomp]
                  confusion.mat.select = get.confusion_matrix(truth = Y[test], predicted = prediction.select)
                  BER.select=get.BER(confusion.mat.select) #Guardo BER de cada n-fold iteracion
                  compMod=ncomp
                  
                  #Cada vez que comience un nuevo modelo (VIP=inicial) me elije las variables seleccionadas en ese
                  #mismo modelo y no en el anterior. De esta forma si pasa en el siguiente
                  #paso por el break, se quedan guardadas las variables con el VIPini de este modelo
                  if (threshold == VIPinicial){
                    select.Var1=colnames(Xsel) #Variables seleccionadas
                    Nvar.select=length(select.Var1) #Numero de variables seleccionadas
                    compModround=round(mean(compMod)) #Componentes del modelo
                    BER.r=mean(BER.select) #BER modelo CON reduccion de variables
                    VIPfinal=threshold
                  }
                  
                  if (BER.select>BER){
                    perform=1
                  }else{
                    perform=0
                  }
                  continuar=c(continuar+perform)
                  
                  if (continuar==oportunidad+1){break} else{
                    
                    if (BER.select <= BER){
                      select.Var1=colnames(Xsel) #Variables seleccionadas
                      Nvar.select=length(select.Var1) #Numero de variables seleccionadas
                      compModround=round(mean(compMod)) #Componentes del modelo
                      BER.r=mean(BER.select) #BER modelo CON reduccion de variables
                      VIPfinal=threshold
                    }
                  }
                  #########################################################################################
                  #Seleccion de VIP del modelo nuevo
                  #########################################################################################
                  my.plsda.perf.select=mixOmics::plsda(Xsel, Ysel, ncomp = round(mean(compMod))) #Para obtener VIP hago modelo con 
                  #la reduccion de variables pero todas las muestras (de manera analoga a antes de la reduccion de variables)
                  my_vip=vip(my.plsda.perf.select)
                  selectVIP=c()
                  
                  for(selC in 1:dim(my_vip)[1]){
                    selectVIP[selC]=my_vip[selC,which.max(my_vip[selC,])]
                  }
                  as.matrix(selectVIP)
                  if (is.null(BER.select)==TRUE){
                  }else{
                    if(BER.select<=BER){
                      BER = BER.select
                    }
                  }
                }
                
                var.select[[cv]] <- list()
                
                var.select[[cv]]$VIPfinal = VIPfinal
                var.select[[cv]]$select.Var1 = select.Var1
                var.select[[cv]]$Nvar.select = Nvar.select
                var.select[[cv]]$compModround = compModround
                var.select[[cv]]$BER.r = BER.r
              }
              
              # parallel::stopCluster(cl)
            }

            message("Obtención de iter.select")
            
            # Obtención de iter.select
            iter.select = iterSelect(X = X, Y = Y, var.select = var.select, 
                                     compPLSDA = compPLSDA, seq.iter = seq.iter, repetcv = repetcv)
            
            #########################################
            # Obtener matriz con las nuevas variables
            #frequenceMinStability = 0.92 #Selecciono las condiciones definitivas
            #########################################
            
            minStability = frequenceMinStability * repetcv
            remainVar = c()
            variables = iter.select$variables
            choise=rownames(variables)[which(variables>minStability)]
            for (choisei in 1:length(choise)){
              remainVar[choisei] = which(colnames(X)==choise[choisei])
            }
            Xnew = X[,remainVar]
            
            # Return ####
            features_selected = colnames(Xnew)
            
            dfData <- object@features %>%
              dplyr::select(all_of(features_selected[order(features_selected)]))
            
            object@features <- dfData %>% as_tibble()
            
            object@metadata$featuresSelectionMethod <- paste("vip ",frequenceMinStability)
            
            # CSV Export #### 
            if(export) {
              ASWExportCSV(object, "output/fs/", "rfe")  
            }   
            
            return (object)
          }
          
)




# Method: ASWVIPbasedFeatureSelectionStrict ####
setGeneric("ASWVIPbasedFeatureSelectionStrict", function(object, verbose=F, export=F, plot=F) {
  standardGeneric("ASWVIPbasedFeatureSelectionStrict")
})

setMethod("ASWVIPbasedFeatureSelectionStrict", "ASWExperiment",
          function(object, verbose, export, plot){
            return (ASWVIPbasedFeatureSelection(object,frequenceMinStability=0.9, verbose, export, plot))
          })
            
# Method: ASWVIPbasedFeatureSelectionRelaxed ####
setGeneric("ASWVIPbasedFeatureSelectionRelaxed", function(object, verbose=F, export=F, plot=F) {
  standardGeneric("ASWVIPbasedFeatureSelectionRelaxed")
})

setMethod("ASWVIPbasedFeatureSelectionRelaxed", "ASWExperiment",
          function(object, verbose, export, plot){
            return (ASWVIPbasedFeatureSelection(object,frequenceMinStability=0.7, verbose, export, plot))
          })

# Method: ASWFFisherFeatureSelection####
setGeneric("ASWFFisherFeatureSelection", function(object, prop=0.2, verbose=F, export=F, plot=F) {
  standardGeneric("ASWFFisherFeatureSelection")
})

setMethod("ASWFFisherFeatureSelection", "ASWExperiment",
          function(object, prop, verbose, export, plot){
            
            if (verbose) {
              cat(primeColor("\n[FEATURES SELECTION] "))
              cat(magenta("Fisher's features selection "))
              cat(red("[prop=",prop,"%]\n", sep = ""))
            }
            
            ##
            
            X <- object@features
            Y <- object@style
            
            #
            getFValue <- function(x) {
              return (anova(lm(x ~ Y))$`F value`[1])
            }
            
            nFeatures <- length(X)
            nFeaturesSelected <- round(nFeatures*prop)
            

            # Calculate F Fisher ####
            ffisher <- list()
            
            for(i in 1:nFeatures){
              ffisher[i] <- getFValue(as.matrix(X[i]))
            }
            
            # Get DF with FValue
            df <- X %>%
              t() %>%
              as.data.frame() %>%
              rowwise() %>%
              mutate(fValue = getFValue(c_across())) 
            
            # Get max values by Feature index
            lst <- sort(df$fValue, index.return=TRUE, decreasing=TRUE)
            lstIndex <- lapply(lst, `[`, lst$x %in% head(unique(lst$x),nFeaturesSelected))  
            
            # Return ####
            selectedFeatures <- sort(names(X[lstIndex$ix]))

            object@features <- X[,selectedFeatures]
            
            object@metadata$featuresSelectionMethod <- "ffisher"
            

            # CSV Export #### 
            if(export) {
              ASWExportCSV(object, "output/fs/", "ffisher")  
            }   
            
            return (object)
          }
          
)



