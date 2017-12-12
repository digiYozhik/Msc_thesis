#'Function to make inference on a cross validation analysis and a multi-location 
#'trial data set using BGLR
#' 
#'@param P a data frame that holds the design information and phenotypes for the 
#'data to be modeled. The data should hold following features, which are 
#'detailed below and represent the columns in the data frame. For inference on 
#'the thesis data set we use the data set P obtained by typing data(P) in the 
#'console.
#'@param id character describing the column name for the names of the observations. 
#'Default is GERMPLASM.
#'@param factor character describing the column name for the factor that describes 
#'the geographic location in the mutli-location trial. Default is LOCATION.
#'@param trait character describing the column name for the phenotype to be 
#'modeled. Default is YIELD.
#'@param CVscheme data frame output from the crossValidate function which was based 
#'on a user decided sampling strategy to use in the cross-validation. Default is NULL, 
#'which applies prediction on full data set, and which is not yet implemented, making 
#'specificiation of this argument required.
#'@param modelName character name describing the model used for modeling:
#'\describe{
#'  \item{\code{modelG}:}{Model where entries and locations are seen as random 
#'  terms in the model. The G-matrix is used to include the genetic relatedness 
#'  between the entries. See reference 4 for more detail.}
#'  \item{\code{modelGE}:}{Model where entries and locations are seen as random 
#'  terms in the model, and where a GxE interaction term is included. The G-matrix 
#'  is used to include the genetic relatedness between the entries. The GxE 
#'  interactions are modeled following Jarquin et al. (2014). See reference 3 and 
#'  4 for more detail.}
#'  \item{\code{modelL}:}{Model where entries and locations are seen as random 
#'  terms in the model, and where no information about the relatedness between 
#'  the entries in included. This model is included for didactic and testing 
#'  purposes.}
#'}  
#'@param G matrix containing the realized G-matrix obtained for the entries in 
#'the dataset specified in P.
#'@param outputDir character specifying the name of the directory where to output 
#'the files used in the modeling and inference. Default is the working directory.
#'@param verbose logical whether to output information about the progress of the 
#'cross-validation. Default is FALSE.
#'@param replications numeric defining the number of replications of the 
#'cross-validation. Default is 3. 
#'@param \dots additional arguments for the BGLR function. Of interest are 
#'nIter for the number of iterations and burnIn specifying the 
#'burn-in used in MCMC analysis.
#'@details The function uses the cross-validation scheme information (\code{CVscheme} 
#'argument) to split the data into training and test sets. While running through 
#'the replications (\code{replication} argument) and folds, the model specified in 
#'the \code{modelName} argument is fitted using the BGLR framework following the 
#'specifications in Appendix B of reference 4. After model fit a series of 
#'metrics are calculated to support inference, which is further detailed in 
#'reference 4. This includes the predictive ability, the mean squared prediction 
#'error (MSPE), and the bias which is calculated using a linear model 
#'(\code{lm} function) of observed phenotype values on the predicted phenotype 
#'values of the test set under evaluation. The outputted information is detailed 
#'in the \code{Value} section. The files used for inference are stored in a folder 
#'named BGLR which is a subdirectory of the directory specified in the 
#'\code{outputDir} argument.
#'@return list with following slots, where TS stands for test set.
#'\itemize{
#'  \item{\code{n.SNP} }{Number of SNPs used in analysis. 
#'  Not relevant here, put to zero.}
#'  \item{\code{n.T} }{Matrix with number of entries in the test set for each fold 
#'  (rows) by replications (columns).}
#'  \item{\code{n.DS} }{Matrix with the number of observations in the total 
#'  dataset for each fold(rows) by replications (columns).}
#'  \item{\code{id.TS} }{List of IDs of each test set within a list of each 
#'  replication.}
#'  \item{\code{bu} }{Estimated fixed and random effects of each fold within 
#'  each replication (see crossVal function)}
#'  \item{\code{y.TS} }{Predicted values of all test sets within each replication.}
#'  \item{\code{PredAbi} }{Predictive ability of each fold within each replication 
#'  calculated as correlation coefficient \eqn{r(y_{TS},\hat y_{TS})}.}
#'  \item{\code{rankCor} }{Spearman's rank correlation of each fold within each 
#'  replication calculated between \eqn{y_{TS}} and \eqn{\hat y_{TS}}.}
#'  \item{\code{bias} }{Regression coefficients of a regression of the observed 
#'  values on the predicted values in the TS. A regression coefficient \eqn{< 1} 
#'  implies inflation of predicted values, and a coefficient of \eqn{> 1} 
#'  deflation of predicted values.}
#'  \item{\code{k} }{Integer defining the number of folds.}
#'  \item{\code{Rep} }{Numeric defining the number of replications.}
#'  \item{\code{sampling} }{Character defining the sampling method.}
#'  \item{\code{Seed} }{Seed for \code{set.seed()}}
#'  \item{\code{rep.seed} }{vector with the values for the seeds used for each 
#'  replication}
#'  \item{\code{nr.ranEff} }{Number of random effects used (see crossVal function)}
#'  \item{\code{VC.est.method} }{Method for the variance components 
#'  (\code{committed} or \code{re-estimated with ASReml/BRR/BL}), 
#'  see crossVal function. We recommend the default, BGLR.}
#'  \item{\code{m10} }{Mean of observed values for the 10\% best predicted of 
#'  each replication. The k test sets are pooled within each replication.}
#'  \item{\code{mse} }{Mean squared error (of prediction, MSPE) of each fold 
#'  within each replication calculated between \eqn{y_{TS}} and \eqn{\hat y_{TS}}. 
#'  This is in reference 4 referred to as MSPE, the mean squared prediction 
#'  error.}
#'  \item{\code{topRecovery} }{Array of topx recovery of entries across the 
#'  different locations. Array contains a matrix for every fold in the cross-
#'  validation. Every matrix hold as as many rows as replications defined. The 
#'  columns in the matrix hold values for the different topx recoveries, where 
#'  x is element of (10, 20, 30, 40, 50, 100, 200). The elements in the matrix 
#'  are calculated as the percentage entries intersecting between the entries in 
#'  the raw and predicted test set under consideration.}
#'  \item{\code{residualErrors} }{Matrix of residual errors, with as columns the 
#'  different folds in the cross-validation and the number of columns representing 
#'  the different replications. The variance was taken from the varE component in 
#'  the fitted BGLR object.}
#'  }
#'@references \describe{
#'      \item{\code{1}:}{Albrecht, T., et al. (2011). Genome-based prediction of 
#'      testcross values in maize. Theor Appl Genet 123:339-350.}
#'      \item{\code{2}:}{De Los Campos, G., Perez, P. (2014). BGLR: Bayesian 
#'      Generalized Linear Regression. Version 1.0.3. 
#'      (http://CRAN.R-project.org/package=BGLR).}
#'      \item{\code{3}:} {Jarquin, D. et al. (2014). A reaction norm model for 
#'      genomic selection using high-dimensional genomic and environmental data. 
#'      Theor Appl Genet 127(3):595-607.}
#'      \item{\code{4}:} {Derijcker, R. (2015). Investigating 
#'      incorporation of genotype x environment interaction (G x E) for genomic 
#'      selection in a practical setting. Unpublished M.Sc. thesis. 
#'      University of Ghent:Belgium.}
#'}
#'@export
#'@author Ruud Derijcker
#'@examples
#'data(G)
#'data(P)
#'scheme <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION", k=5, 
#'                         replication=2, seed=NULL, exclusive=TRUE, 
#'                         sampling="randomByID",verbose=TRUE)
#'output <- inferenceBGLR(P, CVscheme=scheme, modelName="modelG", id="GERMPLASM",
#'                        G=G, factor="LOCATION", trait="YIELD", nIter=1500, burnIn=250,
#'                        replications=2)
#'str(output)
#'
inferenceBGLR <- function(P, id="GERMPLASM", factor="LOCATION", trait="YIELD",
                          CVscheme=NULL, modelName=c("model1","model2", "model0"), 
                          G=G, outputDir=getwd(), verbose=TRUE, replications=3, ...) { 
  if(is.null(CVscheme))
    stop("Full prediction not yet implemented.")
  # prepare data
  if(!is.null(factor)) {
    data <- P[,which(colnames(P) %in% c(id, factor, trait))]
    colnames(data) <- c("IDUnique","FACTOR", "TRAIT")
    data$ID <- paste(data$IDUnique, data$FACTOR, sep="_")
  } else {
    data <- P[,which(colnames(P) %in% c(ID, trait))]
    colnames(data) <- c("ID","TRAIT")
  }
  data <- na.omit(data)
  
  # make folders to store
  setwd(outputDir)
  if (.Platform$OS.type == "unix") {
    BRRTest <- system(paste("ls"), intern = TRUE)
    if (!any(BRRTest %in% "BGLR")) 
      system(paste("mkdir BGLR"))
  }
  if (.Platform$OS.type == "windows") {
    BRRTest <- shell(paste("dir /b"), intern = TRUE)
    if (!any(BRRTest %in% "BGLR")) 
      shell(paste("md BGLR"))
  }
  
  # initiate analysis
  # replications <- ncol(CVscheme)-1
  k <- max(CVscheme[2])
  seed2 <- round(runif(replications, 1, 1e+05), 0)
  topRecovery <- c(10, 20, 30, 40, 50, 100, 200)
  nFactor <- length(unique(P[,colnames(P)[which(colnames(P)==factor)]]))
  foldRepMatrix <- errorMatrix <- matrix(NA, nrow = k, ncol = replications)
  repMatrix <- matrix(NA, nrow=replications, ncol=1)
  predictedMatrix <- matrix(NA, nrow=nrow(P), ncol=replications)
  rownames(predictedMatrix) <- data$ID
  rownames(foldRepMatrix) <- rownames(errorMatrix) <- paste("fold", 1:k, sep = "")
  colnames(foldRepMatrix) <- rownames(repMatrix) <- colnames(errorMatrix) <- colnames(predictedMatrix) <- paste(
    "rep", 1: replications, sep = "")
  correlationPearson <- correlationSpearman <- nTestset <- nDataset <- foldRepMatrix
  bias <- MSE <- m10 <- t(repMatrix)
  t10 <- array(matrix(NA),c(replications,length(topRecovery),nFactor), 
               dimnames=list(c(paste("rep", 1:replications, sep = "")),
                             c(paste("top",topRecovery, sep=""))))
  idTestsetRep <- list()
  predictedByFoldMatrixRep <- NULL
  
  # perform analysis
  for (iReplication in 1:replications) {
    idTestset <- list()
    predictedPhenotype <- rawPhenotype <- NULL
    predictedByFoldMatrix <- matrix(NA, ncol = k, nrow = nrow(data))
    rownames(predictedByFoldMatrix) <- data$ID
    colnames(predictedByFoldMatrix) <- paste("rep", iReplication, 
                                             "_fold", 1:k, sep = "")
    
    for (iFold in 1:k) {
      if (verbose) 
        cat("Replication: ", iReplication, "\t Fold: ", iFold, " \n")
      trainingSet <- CVscheme[!(CVscheme[, 1 + iReplication] %in% iFold), ]
      testSet <- CVscheme[!is.na(CVscheme[, 1 + iReplication]), ]
      testSet <- testSet[testSet[, 1+iReplication] == iFold, ]
      namesDataset <- c(as.character(trainingSet[,"ID"]), 
                        as.character(testSet[,"ID"]))
      dataMaskTestdata <- data
      iskFold <- CVscheme[, 1+iReplication] == iFold
      iskFold[is.na(iskFold)] <- FALSE
      dataMaskTestdata[iskFold == TRUE, "TRAIT"] <- NA
      y <- dataMaskTestdata$TRAIT
      
      if(modelName == "modelG") {
        ETA <- list()
        ETA[[1]] <- list(~factor(dataMaskTestdata$FACTOR)-1, model='BRR')
        L <- t(chol(G))
        indexIDs <- as.integer(factor(x=as.character(dataMaskTestdata$IDUnique),
                                      levels=rownames(G), ordered=TRUE))
        ZL <- L[indexIDs,]
        ETA[[2]]<-list(X=ZL, model='BRR')
        
      } else if (modelName == "modelGE") {
        ETA <- list()
        ZE <- as.matrix(model.matrix(~as.factor(
          dataMaskTestdata$FACTOR)-1)) # Design matrix for environment
        g <- factor(dataMaskTestdata$IDUnique,levels=rownames(G),
                    ordered=TRUE) # factor for germplasm ID
        Zg <- model.matrix(~g-1)
        ZGZ <- Zg %*% G %*% t(Zg)
        EVD.ZGZ <- eigen(ZGZ)
        ExG <- ZGZ*tcrossprod(ZE)
        EVD.ExG<-eigen(ExG)
        ETA[[1]] <- list(X=ZE,model="BRR")
        ETA[[2]]<-list(V=EVD.ZGZ$vectors,d=EVD.ZGZ$values,model='RKHS')
        ETA[[3]]<-list(V=EVD.ExG$vectors,d=EVD.ExG$values,model='RKHS')
        
      } else if (modelName == "modelL") {
        ETA <- list()
        ETA[[1]] <- list(~factor(dataMaskTestdata$FACTOR)-1,model='BRR')
        ETA[[2]]<-list(~factor(dataMaskTestdata$IDUnique)-1,model='BRR')
      }
      capture.output(model <- BGLR(y=y, ETA=ETA,
                                  saveAt = paste("BGLR/rep", iReplication, 
                                                "_fold", iFold, sep = ""), ...), 
                     file = paste("BGLR/BGLRout_rep", iReplication, 
                                  "_fold", iFold, ".txt", sep = ""))
      errorMatrix[iFold, iReplication] <- model$varE
      predictedByFoldMatrix[ ,iFold] <- model$yHat
      yRawTraining <- data[(data$ID %in% 
                              as.character(testSet[, 1])), "TRAIT"] # original
      yPredictedtraining <- model$yHat[(data$ID %in% 
                              as.character(testSet[, 1]))] # predicted
      names(yPredictedtraining) <- names(yRawTraining) <- data$ID[(
        data$ID %in% as.character(testSet[, 1]))]
      nTestset[iFold, iReplication] <- length(yPredictedtraining)
      nDataset[iFold, iReplication] <- length(namesDataset)
      correlationPearson[iFold, iReplication] <- round(cor(yRawTraining,
                       yPredictedtraining), digits = 4)
      correlationSpearman[iFold, iReplication] <- round(cor(yRawTraining, 
                       yPredictedtraining, method = "spearman"), digits = 4)
      predictedPhenotype <- rbind(predictedPhenotype, data.frame(yPredictedtraining))
      rawPhenotype <- rbind(rawPhenotype, data.frame(yRawTraining))
      idTestset[[iFold]] <- as.character(unique(testSet[, 1]))
      names(idTestset)[[iFold]] <- paste("fold", iFold, sep = "")
    }
    
    fitLinearModel <- lm(as.numeric(rawPhenotype$yRawTraining) ~ 
                           as.numeric(predictedPhenotype$yPredictedtraining))
    bias[,iReplication]  <- fitLinearModel$coefficients[2]      
    MSE[,iReplication] <- mean((rawPhenotype - as.numeric(
      predictedPhenotype$yPredictedtraining))^2)
    predictedMatrix[,iReplication] <- predictedPhenotype[order(rownames(
      predictedPhenotype)), ]
    predictedByFoldMatrixRep <- cbind(predictedByFoldMatrixRep, 
                                      predictedByFoldMatrix)
    idTestsetRep[[iReplication]] <- idTestset
    names(idTestsetRep)[[iReplication]] <- paste("rep", iReplication, sep = "")
    
    top10Predicted <- predictedPhenotype$yPredictedtraining[
      order(-predictedPhenotype$yPredictedtraining)]
    names(top10Predicted) <- data$ID
    n10 <- round(0.1 * length(predictedPhenotype))
    top10Predictedsel <- top10Predicted[1:n10]
    m10[,iReplication] <- mean(data[data$ID %in% names(top10Predictedsel), 
                                    "TRAIT"])
    
    if(factor == "LOCATION") {
      df <- data.frame(raw=rawPhenotype, predicted=predictedPhenotype)
      dfNames <- data.frame(unlist(lapply(strsplit(rownames(df), split="_"), 
                                          function(x) x[1])),
                            unlist(lapply(strsplit(rownames(df), split="_"), 
                                          function(x) x[2])))
      colnames(dfNames) <- c("ID","LOCATION")
      df <- data.frame(df, dfNames)
    }
    
    dfList <- dataframe.to.list(df, by="LOCATION")
    rankObserved <- lapply(dfList, function(x) 
      as.character(x$ID[order(x$yRawTraining, decreasing=TRUE)]))
    rankPredicted <- lapply(dfList, function(x) 
      as.character(x$ID[order(x$yPredictedtraining, decreasing=TRUE)]))
    
    for(mRecovery in 1:length(topRecovery)){
      idx <- topRecovery[mRecovery]
      for(n in 1: length(rankObserved)){
        t10[iReplication, mRecovery, n] <- length(intersect(
          rankObserved[[n]][1:idx], rankPredicted[[n]][1:idx]))/idx
      }
    }
  }
  
  obj <- list(n.SNP=0, n.T =nTestset, n.DS=nDataset, 
              id.TS=idTestsetRep, bu=predictedByFoldMatrixRep, 
              y.TS=predictedMatrix, PredAbi=correlationPearson, 
              rankCor=correlationSpearman, bias=bias, k=k, Rep=replications, 
              sampling="", Seed=NULL, rep.seed=seed2, nr.ranEff = 0, 
              VC.est.method="BGLR", m10=m10, mse=MSE, topRecovery=t10, 
              residualErrors=errorMatrix)
  return(obj)
}

#' Summarize the inference results from  from inferenceBGLR
#' 
#' Summarize the inference results from  from inferenceBGLR
#' 
#'@param P a data frame that holds the design information and phenotypes for the 
#'data to be modeled. The data should hold following features, which are 
#'detailed below and represent the columns in the data frame. For inference on 
#'the thesis data set we use the data set P obtained by typing data(P) in the 
#'console.
#'@param obj object obtained from the inferenceBGLR function.
#'@return A list of summarized inference measures from \code{obj}.
#'@export
#'@author Ruud Derijcker
#'@details The function summarizes the obtained results for the bias, MSPE, and 
#'predictive ability ("correlation") using the average across the different folds 
#'and replications as  was specified in the cross-validation and \code{inferenceBGLR} 
#'function. Confidence intervals for the predictive ability were calculated using 
#'rho +/- 1.96*sqrt((1-rho^2)/(nrow(P)-2)). The confidence intervals for the bias 
#'and MSPE were calculated using a minimum and maximum value calculated as 
#'mean(metric) +/- SE(metric), following \code{synbreed} reporting by T. Albrecht 
#'in the \code{summary.cvData} function. See also reference.
#'@references Wimmer, V., et al. (2012) synbreed: a framework for the analysis 
#'of genomic prediction data using R. Bioinformatics, 28: 2086-2087.
#'@examples 
#'data(G)
#'data(P)
#'scheme <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION", k=5, 
#'                         replication=2, seed=NULL, exclusive=TRUE, 
#'                         sampling="randomByID",verbose=TRUE)
#'output <- inferenceBGLR(P, CVscheme=scheme, modelName="modelG", id="GERMPLASM",
#'                        G=G, factor="LOCATION", trait="YIELD", nIter=1500, burnIn=250,
#'                        replications=2)
#'output <- summarizeInference(P, output)
#'output
#' 
summarizeInference <- function (P, obj) {
  ans <- list()
  ans$k <- obj$k
  ans$Rep <- obj$Rep
  ans$testSets <- obj$n.T
  colmean.b <- colMeans(obj$bias)
  colmean.mse <- colMeans(obj$mse)
  se.b <- sd(colmean.b)/sqrt(length(colmean.b))
  ans$bias <- paste(round(mean(obj$bias), 3), "CI[",round(mean(obj$bias) - se.b, 3),
                    ";",round(mean(obj$bias) + se.b, 3), "]")
  
  se.mse <- sd(colmean.mse)/sqrt(length(colmean.mse))
  ans$mse <- paste(round(mean(sqrt(obj$mse)), 3), "CI[",
                   round(mean(sqrt(obj$mse)) - sqrt(se.mse), 3),
                   ";",round(mean(sqrt(obj$mse)) + sqrt(se.mse), 3), "]")
  ans$Seed <- obj$Seed
  ans$rep.seed <- obj$rep.seed
  ans$topRecovery <- format(apply(obj$topRecovery,2, mean), digits=3, nsmall=3)
  tmp <- obj$PredAbi
  tmp <- cbind(tmp,rowMeans(tmp,na.rm=TRUE))
  tmp <- rbind(tmp,colMeans(tmp,na.rm=TRUE))
  rho  <- tmp[dim(obj$topRecovery)[3],ans$Rep+1]
  li <- rho - 1.96*sqrt((1-rho^2)/(nrow(P)-2))
  ls <- rho + 1.96*sqrt((1-rho^2)/(nrow(P)-2))
  ans$PredAbi <- paste("CV", round(rho,3), " [",round(li,3),";",round(ls,3),"]")
  return(ans)
}
