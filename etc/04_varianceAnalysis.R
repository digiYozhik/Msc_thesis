#===============================================================================
# Analysis script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)             #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#                                                                              #
# Step 4 : - Variance analysis on full data set                                #
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
# 1 - Load and prepare data
#-------------------------------------------------------------------------------
data(M) # genotypes
data(P) # phenotypes
data(G) # realization of the G-matrix

# Load provided object or stored objects from previous scripts
load(file=file.path(inputPath, "ZE.rda"))
load(file=file.path(inputPath, "Zg.rda"))
load(file=file.path(inputPath, "ZGZ.rda"))
load(file=file.path(inputPath, "EVD.ZGZ.rda"))
load(file=file.path(inputPath, "ExG.rda"))
load(file=file.path(inputPath, "EVD.ExG.rda"))
load(file=file.path(inputPath, "L.rda"))
load(file=file.path(inputPath, "MCMCparameters.rda"))

# Scale the yield performance measures prior to analysis
P$YIELDScaled <- scale(P$YIELD,center=TRUE,scale=TRUE)

#-------------------------------------------------------------------------------
# 2 - Variance analysis for models excluding genetic relationships
#-------------------------------------------------------------------------------
# Model with only location term
ETA <- list()
ETA[[1]] <- list(~factor(P$LOCATION)-1,model='BRR')
fmE <- BGLR(y=P$YIELDScaled,ETA=ETA)

# Model with location and entry term
ETA[[2]]<-list(~factor(P$GERMPLASM)-1,model='BRR')
fmEV <- BGLR(y=P$YIELDScaled,ETA=ETA)

#-------------------------------------------------------------------------------
# 3 - Variance analysis for modelG
#-------------------------------------------------------------------------------
L <- t(chol(G))
indexIDs <- as.integer(factor(x=as.character(P$GERMPLASM),levels=rownames(G),
                            ordered=TRUE))
ZL <- L[indexIDs,]
ETA[[2]]<-list(X=ZL,model='BRR')
fmG <- BGLR(y=P$YIELDScaled,ETA=ETA)

#-------------------------------------------------------------------------------
# 4 - Variance analysis for model GE
#-------------------------------------------------------------------------------
ETA[[1]] <- list(X=ZE,model="BRR")
ETA[[2]]<-list(V=EVD.ZGZ$vectors,d=EVD.ZGZ$values,model='RKHS')
ETA[[3]]<-list(V=EVD.ExG$vectors,d=EVD.ExG$values,model='RKHS')
fmGE <- BGLR(y=P$YIELDScaled, ETA=ETA)

#-------------------------------------------------------------------------------
# 5 - Summarize variance components
#-------------------------------------------------------------------------------
resultsVarianceAnalysis <- data.frame(
           varE=c(fmEV$ETA[[1]]$varB, fmG$ETA[[1]]$varB, fmGE$ETA[[1]]$varB),
           varL=c(fmEV$ETA[[2]]$varB, NA, NA),
           varG=c(NA, fmG$ETA[[2]]$varB, fmGE$ETA[[2]]$varU),
           varGxE=c(NA, NA, fmGE$ETA[[3]]$varU),
           varRes=c(fmEV$varE, fmG$varE, fmGE$varE))
rownames(resultsVarianceAnalysis) <- c("modelEV", "modelG", "modelGE")
round(resultsVarianceAnalysis, 2)

#-------------------------------------------------------------------------------
# 6 - Summarize model fitting
#-------------------------------------------------------------------------------
fmEV$fit$DIC
fmG$fit$DIC
fmGE$fit$DIC

#-------------------------------------------------------------------------------
# 7 - Heritability calculations using residual variance
#-------------------------------------------------------------------------------
varE <- resultsVarianceAnalysis$varE
names(varE)<-c('LOCATION','LOCATION + Var(I)','LOCATION + Var(G)')
h2 <- c(fmEV$ETA[[2]]$varB/(fmEV$ETA[[2]]$varB+fmEV$varE), 
      mean(diag(G))*fmG$ETA[[2]]$varB/(mean(diag(G))*fmG$ETA[[2]]$varB+fmG$varE))
names(h2) <- c('LOCATION + Var(I)','LOCATION + Var(G)'); h2
barplot(h2,col=8,ylab='Heritability')

#===============================================================================
# End of script                                                                #
#===============================================================================