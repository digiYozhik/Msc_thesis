#===============================================================================
# Analysis script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)             #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#                                                                              #
# Step 3 : - Align phenotypes and genotypes                                    #
#          - Prepare matrices for modeling                                     #
#          - Preliminary modeling                                              #
#          - MCMC diagnostics                                                  #
#          - Assess prior sensitivity for hyper-prior of variances             #
#===============================================================================
#-------------------------------------------------------------------------------
# DISCLAIMER                                                                   #
# The analysis scripts are used in combination with the user-defined package   #
# MaStatThesisRuud, which is used to obtain results and inference in the M.Sc. # 
# thesis. The package uses the packages BGLR and synbreed as dependencies, and # 
# comes with fully documented functions with (mostly) included examples.       #
# Contact: Ruud Derijcker <ruud.derijcker@ugent.be>                            #
# Some notes about the provided analysis scripts:                              #
# - Scripts are preferentially executed using the ordering by their number     #
#   prefix. Though, scripts can be executed in a stand-alone fashion. For this #
#   purpose, we provide the user with the adequate R-objects                   #
# - Scripts used in preparation of the raw data are not provided because of    # 
#   confidentiality reasons                                                    #
# - Scripts involving phenotypic modeling are shown assuming access to the     #
#   software ASReml-R.                                                         #
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 0 - Prepare data analysis
#-------------------------------------------------------------------------------
rm(list=ls()) # clear all R-objects in memory
# Data directories are organized using a main project folder, with subfolders 
# for the input data and output data:
# projectPath
#    |___ inputPath
#    |___ outputPath
# The user needs to update these folder paths prior to analysis
projectPath <- file.path(getwd(), "thesisProject")
inputPath <- file.path(projectPath, "input")
outputPath <- file.path(projectPath, "output")

# Load needed libraries
library(MaStatThesisRuud) # base funtionalities
help("MaStatThesisRuud")
help(package="MaStatThesisRuud")

#-------------------------------------------------------------------------------
# 1 - Load prepared data
#-------------------------------------------------------------------------------
data(Mraw) # genotypes
data(Praw) # phenotypes
data(Graw) # G-matrix

#-------------------------------------------------------------------------------
# 2 - Align genotypes and phenotypes
#-------------------------------------------------------------------------------
# Remove genotypes not connected to phenotypes
commonNames <- intersect(colnames(Mraw), Praw$GERMPLASM)
(length(commonNames)) # number of unique plants
M <- Mraw[,colnames(Mraw) %in% commonNames]
P <- droplevels(Praw[Praw$GERMPLASM %in% commonNames,])
G <- Graw[rownames(Graw) %in% commonNames, colnames(Graw) %in% commonNames]
dim(M); dim(P); dim(G)

save(M, file=file.path(dataPath, "M.rda"))
save(P, file=file.path(dataPath, "P.rda"))
save(G, file=file.path(dataPath, "G.rda"))

# Alternatively, one can use
data(M) # aligned genotypes
data(P) # aligned phenotypes
data(G) # final G-matrix

#-------------------------------------------------------------------------------
# 3 - Prepare matrices
#-------------------------------------------------------------------------------
#******************************************************************************#
# Prior to modeling, we prepare the matrices needed to do modeling using BGLR. #
# For more detail on the matrices, we refere to Appendix B in the thesis.      #
# These objects are also provided as R-objects (*.rda) to the user             #
#******************************************************************************#
# Sort information based on genotype names
index <- order(colnames(G)) 
G <- G[index,index]
index <- order(as.character(P$GERMPLASM),P$LOCATION)
P <- P[index,]
IDs <- as.character(unique(P$GERMPLASM))

# Prepare incidence matrices
## Indices matrix for locations
ZE <- as.matrix(model.matrix(~as.factor(P$LOCATION)-1))

## Indices matrix for entries
g <- factor(P$GERMPLASM,levels=rownames(G),ordered=TRUE) # factor for germplasm ID
Zg <- model.matrix(~g-1)
if(any(colnames(Zg)!=paste("g",colnames(G),sep=""))) 
  stop("Ordering problem\n")

## Capture genetic information
ZGZ <- Zg %*% G %*% t(Zg)
EVD.ZGZ <- eigen(ZGZ) # Eigen value decomposition: can take a while to calculate!

## Capture G x E information
ExG <- ZGZ*tcrossprod(ZE)
EVD.ExG <- eigen(ExG) # Eigen value decomposition: can take a while to calculate!

## Cholesky decomposition for G-matrix
L <- t(chol(G))

# Save the objects
save(ZE, file=file.path(outputPath, "ZE.rda"))
save(Zg, file=file.path(outputPath, "Zg.rda"))
save(ZGZ, file=file.path(outputPath, "ZGZ.rda"))
save(EVD.ZGZ, file=file.path(outputPath, "EVD.ZGZ.rda"))
save(ExG, file=file.path(outputPath, "ExG.rda"))
save(EVD.ExG, file=file.path(outputPath, "EVD.ExG.rda"))
save(L, file=file.path(outputPath, "L.rda"))

#-------------------------------------------------------------------------------
# 4 - Set MCMC parameters
#-------------------------------------------------------------------------------
# We use following defaults for MCMC sampling in our preliminary modeling
MCMCparameters <- list(nIterations=5000, # number of iterations
                       burnIn=1500, # burn-in
                       thinning=2) # thinning factor
save(MCMCparameters, file=file.path(outputPath, "MCMCparameters.rda"))

#-------------------------------------------------------------------------------
# 5 - Do preliminary modeling
#-------------------------------------------------------------------------------
#******************************************************************************#
# For modeling we use the BGLR package, which provides a General Bayesian      #
# framework for genomic modeling. For details we refer to Appendix B in the    # 
# thesis, and to the BGLR vignette (especially for hyper-prior definitions and #
# rules)                                                                       #
#******************************************************************************#
vignette("BGLR-extdoc")
oldPath <- getwd()
preliminaryPath <- file.path(outputPath, "preliminaryModeling") # set output path
if (file.exists(preliminaryPath)){
  setwd(preliminaryPath)
} else {
  dir.create(preliminaryPath)
  setwd(preliminaryPath)
}

# Model building
## Model with only environment
ETA <- list()
ETA[[1]] <- list(~factor(P$LOCATION)-1,model='BRR')
fmE <- BGLR(y=P$YIELD,ETA=ETA,saveAt="modelE_")

## Model with environment and entries (no relatedness information used)
ETA[[2]]<-list(~factor(P$GERMPLASM)-1,model='BRR')
fmEV <- BGLR(y=P$YIELD,ETA=ETA, saveAt="modelEV_")

## Model G (y = \mu + g + e)
L <- t(chol(G))
indexIDs <- as.integer(factor(x=as.character(P$GERMPLASM),levels=rownames(G),
                              ordered=TRUE))
ZL <- L[indexIDs,]
ETA[[2]]<-list(X=ZL,model='BRR') # BRR: Bayesian Ridge Regression
# Next line takes a while to calculate ...
fmEG <- BGLR(y=P$YIELD,ETA=ETA, saveAt="modelG_",thin=1, nIter=40000, 
             burnIn=1000) # We use higher MCMC numbers since use in diagnostics

## Model GE (y = X\beta + G/Z + e with G/Z:genotypes + gxE)
ETA[[1]] <- list(X=ZE,model="BRR")
ETA[[2]]<-list(V=EVD.ZGZ$vectors, d=EVD.ZGZ$values, model='RKHS')
ETA[[3]]<-list(V=EVD.ExG$vectors, d=EVD.ExG$values, model='RKHS')
fmExG <- BGLR(y=P$YIELD, ETA=ETA, saveAt="modelGE_")

#-------------------------------------------------------------------------------
# 6 - MCMC Diagnostics
#-------------------------------------------------------------------------------
# Read-in data for doing MCMC diagnostics
setwd(preliminaryPath)
fileNames <- list.files(pattern="modelG_")
chains <- matrix(nrow=40000, # MCMCparameters$nIterations-MCMCparameters$burnIn, 
                 ncol=length(fileNames))
colnames(chains) <- fileNames

for(files in 1:length(fileNames)){
  tmp <- read.table(fileNames[files], header=FALSE, sep=" ")
  chains[,files] <- tmp[,1]
}
setwd(oldPath)

# MCMC diagnostics
for(iChain in 1:ncol(chains)){
  tmp <- as.mcmc(chains[ ,iChain])
  plotName <- gsub(".dat","", colnames(chains)[iChain])
  # Trace plot
  traceplot(tmp, smooth=TRUE, main=paste("Traceplot\n",plotName))
  abline(h=mean(chains[ ,iChain]), col="green")
  # Autocorrelation plot
  autocorr.plot(tmp, main=paste("Autocorrelation\n",plotName))
  # Raftery and Lewis's diagnostic
  raft <- raftery.diag(tmp, q=0.5, r=0.005, s=0.95)
  print(raft)
  raft$resmatrix
}
# We assume different convergence criteria for different variables to estimate.
# Following Raftery and Lewis's diagnostic, and some empirical checking of the 
# number of chains needed, we decided to go for 40000 iterations

# Raftery and Lewis's diagnostic output for different variables
# Quantile (q) = 0.5
# Accuracy (r) = +/- 0.005
# Probability (s) = 0.95 
# 
# Burn-in  Total  Lower bound  Dependence
# (M)      (N)    (Nmin)       factor (I)
# 10       204110 38415        5.31      
# 
# 
# Quantile (q) = 0.5
# Accuracy (r) = +/- 0.005
# Probability (s) = 0.95 
# 
# Burn-in  Total  Lower bound  Dependence
# (M)      (N)    (Nmin)       factor (I)
# 25       307540 38415        8.01      
# 
# 
# Quantile (q) = 0.5
# Accuracy (r) = +/- 0.005
# Probability (s) = 0.95 
# 
# Burn-in  Total    Lower bound  Dependence
# (M)      (N)      (Nmin)       factor (I)
# 936      11315556 38415        295       
# 
# 
# Quantile (q) = 0.5
# Accuracy (r) = +/- 0.005
# Probability (s) = 0.95 
# 
# Burn-in  Total Lower bound  Dependence
# (M)      (N)   (Nmin)       factor (I)
# 4        81266 38415        2.12      

#-------------------------------------------------------------------------------
# 7 - Assess prior sensitivity for hyper-prior for variances
#-------------------------------------------------------------------------------
oldPath <- getwd()
priorSensitivityPath <- file.path(outputPath, "priorSensitivity") # set output path
if (file.exists(priorSensitivityPath)){
  setwd(priorSensitivityPath)
} else {
  dir.create(priorSensitivityPath)
  setwd(priorSensitivityPath)
}

#_______________________________________________________________________________
# 7.1 - Settings
#_______________________________________________________________________________
# Define the fold changes we want for the testing of the parameter
folds <- c(0.1, 10, 100)

# Define MCMC parameters (in case not yet defined)
MCMCparameters <- list(nIterations=5000, # number of iterations
                       burnIn=1500, # burn-in
                       thinning=2) # thinning factor
#_______________________________________________________________________________
# 7.2 - Set up Null Model (model G)
#_______________________________________________________________________________
# Fit Null Model
ETA <- list()
ETA[[1]] <- list(~factor(P$LOCATION)-1,model='BRR')
L <- t(chol(G))
indexIDs <- as.integer(factor(x=as.character(P$GERMPLASM),levels=rownames(G),
                              ordered=TRUE))
ZL <- L[indexIDs,]
ETA[[2]]<-list(X=ZL,model='BRR')
NullModel <- BGLR(y=P$YIELD, ETA=ETA, nIter=MCMCparameters$nIterations, 
                  burnIn=MCMCparameters$burnIn, thin=MCMCparameters$thinning,
                  saveAt="NullModel_")
# Extract priors
NullPriors <- list(df0=NullModel$ETA[[1]]$df0, S0=NullModel$ETA[[1]]$S0, 
                   R2=NullModel$ETA[[1]]$R2)
#_______________________________________________________________________________
# 7.3 - Run Null Model with different folds for priors
#_______________________________________________________________________________
for(newPrior in folds) {
  priors <- list(df0=NullPriors$df0, S0=NullPriors$S0*newPrior,
                 R2=NullPriors$R2)
  foldModel <- BGLR(y=P$YIELD, ETA=ETA, nIter=MCMCparameters$nIterations, 
                    burnIn=MCMCparameters$burnIn, thin=MCMCparameters$thinning,
                    df0=priors$df0, S0=priors$S0, R2=priors$R2,
                    saveAt=paste("foldModel_", newPrior, "_", sep=''))
}

#_______________________________________________________________________________
# 7.4 - Plot results for posterior distribution & expected distribution for 
# error variance
#_______________________________________________________________________________
# Extract varE from the different models
varE <- read.table(paste("NullModel_", "varE.dat", sep=""), header=FALSE)  
for(newPrior in folds)
  varE <- data.frame(varE, read.table(paste("foldModel_", newPrior, 
                                            "_varE.dat", sep=''), header=FALSE))
varE <- varE[NullModel$burnIn:nrow(varE), ]

# Prepare density plots
DensityForPlot <- list("vector", length=ncol(varE))
for (model in 1:ncol(varE)) { 
  DensityForPlot[[model]] <- density(varE[,model]) 
  if(!exists("XMax")) { 
    XMax <- max(DensityForPlot[[model]]$x)
    YMax <- max(DensityForPlot[[model]]$y)
  } else {
    XMax <- ifelse(max(DensityForPlot[[model]]$x)>XMax, 
                   max(DensityForPlot[[model]]$x), XMax)
    YMax <- ifelse(max(DensityForPlot[[model]]$y)>YMax, 
                   max(DensityForPlot[[model]]$y), YMax)
  }
}
# xMax <- 60000

# Plot posterior distribution
png(file.path(outputPath, "sensitivityE.png"))
  plot(1, type="n", xlim=c(0,XMax), ylim=c(0, YMax),
       main="Prior Sensitivity Error Variance")
  colors <- c("black", "red", "blue", "green", "magenta", "indianred2", 
              "lightgreen")
  for (plots in 1:length(DensityForPlot))
    points(DensityForPlot[[plots]], type="l", col=colors[plots])
  legend(x="topright", legend=c("NullModel", folds), bty="n", col=colors, pch=3)
  
  # Plot expected distribution
  expected <- 1/(rgamma(NullModel$nIter, 
                        shape=(NullPriors$df0/2), 
                        rate=(2*(NullPriors$S0*NullPriors$df0))))
  points(density(expected), type="l", lty="dashed", col=colors[1])
  
  for(expected in 1:length(folds)){
    expected <- 1/(rgamma(NullModel$nIter, 
                      shape=(NullPriors$df0/2), 
                      rate=(2*(NullPriors$S0*NullPriors$df0*folds[expected]))))
    points(density(expected), type="l", lty="dashed", col=colors[expected+1])
  }
dev.off()
rm(varE, XMax, YMax, DensityForPlot)

#_______________________________________________________________________________
# 7.5 - Plot results for posterior distribution & expected distribution for 
# entries variance
#_______________________________________________________________________________
# Extract varM from the different models
varM <- read.table(paste("NullModel_",  "ETA_1_", "varB",".dat", sep=""), 
                   header=FALSE)  
for(newPrior in folds)
  varM <- data.frame(varM, read.table(paste("foldModel_", newPrior, 
                  "_ETA_1_", "varB",".dat", sep=''), header=FALSE))
varM <- varM[NullModel$burnIn:nrow(varM), ]

# Prepare density plots
DensityForPlot <- list("vector", length=ncol(varM))
for (model in 1:ncol(varM)) { 
  DensityForPlot[[model]] <- density(varM[,model]) 
  if(!exists("XMax")) { 
    XMax <- max(DensityForPlot[[model]]$x)
    YMax <- max(DensityForPlot[[model]]$y)
  } else {
    XMax <- ifelse(max(DensityForPlot[[model]]$x)>XMax, 
                   max(DensityForPlot[[model]]$x), XMax)
    YMax <- ifelse(max(DensityForPlot[[model]]$y)>YMax, 
                   max(DensityForPlot[[model]]$y), YMax)
  }
}

# Plot posterior distribution
png(file.path(outputPath, "sensitivityB.png"))
  plot(1, type="n", xlim=c(0,XMax), ylim=c(0, YMax),
       main="Prior Sensitivity Model Variance")
  colors <- c("black", "red", "blue", "green", "magenta", 
              "indianred2", "lightgreen")
  for (plots in 1:length(DensityForPlot))
    points(DensityForPlot[[plots]], type="l", col=colors[plots])
  legend(x="topright", legend=c("NullModel", folds), bty="n", col=colors, pch=3)
  
  ## plot expected distribution
  expected <- 1/(rgamma(NullModel$nIter, 
                        shape=(NullPriors$df0/2), 
                        rate=(2*(NullPriors$S0*NullPriors$df0))))
  points(density(expected), type="l", lty="dashed", col=colors[1])
  for(expected in 1:length(folds)){
    expected <- 1/(rgamma(NullModel$nIter, 
                      shape=(NullPriors$df0/2), 
                      rate=(2*(NullPriors$S0*NullPriors$df0*folds[expected]))))
    points(density(expected), type="l", lty="dashed", col=colors[expected+1])
  }
dev.off()

#===============================================================================
# End of script                                                                #
#===============================================================================