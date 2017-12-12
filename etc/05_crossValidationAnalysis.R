#===============================================================================
# Analysis script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)             #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#                                                                              #
# Step 5 : - Prepare cross-validation schemes                                  #
#          - Do cross-validation analysis                                      #
#          - Summarize inference                                               #
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
# 1 - Load data
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

#-------------------------------------------------------------------------------
# 2 - Illustrate different CV-schemes
#-------------------------------------------------------------------------------
data(exampleCV)
help(crossValidate)
example(plotCV)

# plots for thesis
germplasm <- c(paste("PLANT0", rep(1:9), sep=""), 
               paste("PLANT", rep(10:20), sep=""))
data <- data.frame(GERMPLASM= germplasm, 
                   LOCATION=rep(1:5, each=20), trait=rnorm(100, 0,1),
                   RANGE=rep(1:4, each=5), ROW=rep(1:5, each=1))
data <- data[order(data$GERMPLASM, data$LOCATION),]
data$LOCATION <- as.factor(paste("Location", data$LOCATION, sep=""))

title <- "Field Design 1"
data1 <- data.frame(data, scheme1)
data1 <- data1[order(data1$ROW, data1$RANGE),]
data1$"Phenotyped" <- as.factor(ifelse(data1$Rep1==1,"no", "yes"))
pg1 <- ggplot(data1, aes(ROW, RANGE, fill = Phenotyped)) +
  facet_grid(~ LOCATION) + geom_tile() + labs(title=title) +
  scale_fill_manual(values = c("grey","green3"))
plot(pg1)

title <- "Field Design 2"
data2 <- data.frame(data, scheme2)
data2 <- data2[order(data2$ROW, data2$RANGE),]
data2$"Phenotyped" <- as.factor(ifelse(data2$Rep1==1,"no", "yes"))
pg2 <- ggplot(data2, aes(ROW, RANGE, fill = Phenotyped)) +
  facet_grid(~ LOCATION) + geom_tile() + labs(title=title) +
  scale_fill_manual(values = c("grey","green3"))
plot(pg2)

title <- "Field Design 3"
data3 <- data.frame(data, scheme3)
data3 <- data3[order(data3$ROW, data3$RANGE),]
data3$"Phenotyped" <- as.factor(ifelse(data3$Rep1==1,"no", "yes"))
pg3 <- ggplot(data3, aes(ROW, RANGE, fill = Phenotyped)) +
  facet_grid(~ LOCATION) + geom_tile() + labs(title=title) +
  scale_fill_manual(values = c("grey","green3"))
plot(pg3)

png(file.path(outputPath, "fieldDesigns.png"), width = 1200,height = 400)
grid.arrange(pg1, pg2, pg3, ncol=3)
dev.off()

#-------------------------------------------------------------------------------
# 3 - Prepare CV-schemes on data
#-------------------------------------------------------------------------------
# Prepare the CV-schemes
scheme1 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=10, seed=NULL, exclusive=TRUE, 
                         sampling="randomByFactor",verbose=TRUE)  
scheme2 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=10, seed=NULL, exclusive=TRUE, 
                         sampling="incompleteTrial",verbose=TRUE)
scheme3 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=10, seed=NULL, exclusive=TRUE, 
                         sampling="randomAccrossFactor",verbose=TRUE)  
scheme4 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=10, seed=NULL, exclusive=TRUE, 
                         sampling="randomWithinFactor",verbose=TRUE)
scheme5 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=10, seed=NULL, exclusive=TRUE, 
                         sampling="randomByID",verbose=TRUE)  

# Plot the CV-schemes
png(file.path(texFigurePath, "scheme3.png"))
plotCVThesisData(P, scheme3, "CV scheme 3")
dev.off()

png(file.path(texFigurePath, "scheme4.png"))
plotCVThesisData(P, scheme4, "CV scheme 4")
dev.off()

png(file.path(texFigurePath, "scheme5.png"))
plotCVThesisData(P, scheme5, "CV scheme 5")
dev.off()

png(file.path(texFigurePath, "scheme1.png"))
plotCVThesisData(P, scheme1, "CV scheme 1")
dev.off()

png(file.path(texFigurePath, "scheme2.png"))
plotCVThesisData(P, scheme2, "CV scheme 2")
dev.off()

#-------------------------------------------------------------------------------
# 4 - Cross-validation analysis on modelG
#-------------------------------------------------------------------------------
# This step can take several hours!
GBLUP_CV1 <- inferenceBGLR(P, CVscheme=scheme1, modelName="modelG", id="GERMPLASM", 
                           G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                           burnIn=10000, replications=10)
save(GBLUP_CV1, file=file.path(inputPath, "GBLUP_CV1.rda"))

GBLUP_CV2 <- inferenceBGLR(P, CVscheme=scheme2, modelName="modelG", id="GERMPLASM", 
                           G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                           burnIn=10000, replications=10)
save(GBLUP_CV2, file="GBLUP_CV2.rda")

GBLUP_CV3 <- inferenceBGLR(P, CVscheme=scheme3, modelName="modelG", id="GERMPLASM", 
                           G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                           burnIn=10000, replications=10)
save(GBLUP_CV3, file="GBLUP_CV3.rda")

GBLUP_CV4 <- inferenceBGLR(P, CVscheme=scheme4, modelName="modelG", id="GERMPLASM", 
                           G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                           burnIn=10000, replications=10)
save(GBLUP_CV4, file="GBLUP_CV4.rda")

GBLUP_CV5 <- inferenceBGLR(P, CVscheme=scheme5, modelName="modelG", id="GERMPLASM", 
                            G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                            burnIn=10000, replications=10)
save(GBLUP_CV5, file="GBLUP_CV5.rda")

#-------------------------------------------------------------------------------
# 5 - Summarize inference for modelG
#-------------------------------------------------------------------------------
#******************************************************************************#
# The user can load the result files from the provided objects as shown below  #
#******************************************************************************#
load(file=file.path(inputPath, "GBLUP_CV1.rda"))
load(file=file.path(inputPath, "GBLUP_CV2.rda"))
load(file=file.path(inputPath, "GBLUP_CV3.rda"))
load(file=file.path(inputPath, "GBLUP_CV4.rda"))
load(file=file.path(inputPath, "GBLUP_CV5.rda"))

# Summarize average predictive abilities for the different CV-schemes
mean(GBLUP_CV1$PredAbi)
mean(GBLUP_CV2$PredAbi)
mean(GBLUP_CV3$PredAbi)
mean(GBLUP_CV4$PredAbi)
mean(GBLUP_CV5$PredAbi)

# Summarize inference for the different CV-schemes
summarizeInference(P, GBLUP_CV1)
summarizeInference(P, GBLUP_CV2)
summarizeInference(P, GBLUP_CV3)
summarizeInference(P, GBLUP_CV4)
summarizeInference(P, GBLUP_CV5)

#-------------------------------------------------------------------------------
# 6 - Cross-validation analysis on modelGE
#-------------------------------------------------------------------------------
# This step can take several hours!
RKHS_CV1 <- inferenceBGLR(P, CVscheme=scheme1, modelName="model2", id="GERMPLASM", 
                          G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                          burnIn=10000, replications=10)
save(RKHS_CV1, file="RKHS_CV1.rda")

RKHS_CV2 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model2", id="GERMPLASM", 
                          G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                          burnIn=10000, replications=10)
save(RKHS_CV2, file="RKHS_CV2.rda")

RKHS_CV3 <- inferenceBGLR(P, CVscheme=scheme3, modelName="model2", id="GERMPLASM", 
                          G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                          burnIn=10000, replications=10)
save(RKHS_CV3, file="RKHS_CV3.rda")

RKHS_CV4 <- inferenceBGLR(P, CVscheme=scheme4, modelName="model2", id="GERMPLASM", 
                          G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                          burnIn=10000, replications=10)
save(RKHS_CV4, file="RKHS_CV4.rda")

RKHS_CV5 <- inferenceBGLR(P, CVscheme=scheme5, modelName="model2", id="GERMPLASM", 
                          G=G, factor="LOCATION", trait="YIELD", nIter=40000, 
                          burnIn=10000, replications=10)
save(RKHS_CV5, file="RKHS_CV5.rda")

#-------------------------------------------------------------------------------
# 7 - Summarize inference for modelG
#-------------------------------------------------------------------------------
#******************************************************************************#
# The user can load the result files from the provided objects as shown below  #
#******************************************************************************#
load(file=file.path(inputPath, "RKHS_CV1.rda"))
load(file=file.path(inputPath, "RKHS_CV2.rda"))
load(file=file.path(inputPath, "RKHS_CV3.rda"))
load(file=file.path(inputPath, "RKHS_CV4.rda"))
load(file=file.path(inputPath, "RKHS_CV5.rda"))

# Summarize average predictive abilities for the different CV-schemes
mean(RKHS_CV1$PredAbi)
mean(RKHS_CV2$PredAbi)
mean(RKHS_CV3$PredAbi)
mean(RKHS_CV4$PredAbi)
mean(RKHS_CV5$PredAbi)

# Summarize inference for the different CV-schemes
summarizeInference(P, RKHS_CV1)
summarizeInference(P, RKHS_CV2)
summarizeInference(P, RKHS_CV3)
summarizeInference(P, RKHS_CV4)
summarizeInference(P, RKHS_CV5)

#-------------------------------------------------------------------------------
# 8 - Summarize inference for difference of modelG and modelGE
#-------------------------------------------------------------------------------
# CV1
cat("[", round(mean(RKHS_CV1$PredAbi-GBLUP_CV1$PredAbi)-sd(diff)/sqrt(5), 4), ";", 
    round(mean(RKHS_CV1$PredAbi-GBLUP_CV1$PredAbi) + (sd(diff)/sqrt(5)), 4),"]") 

# CV2
cat("[", round(mean(RKHS_CV2$PredAbi-GBLUP_CV2$PredAbi)-sd(diff)/sqrt(5), 4), ";", 
    round(mean(RKHS_CV2$PredAbi-GBLUP_CV2$PredAbi) + (sd(diff)/sqrt(5)), 4),"]") 

# CV3
cat("[", round(mean(RKHS_CV3$PredAbi-GBLUP_CV3$PredAbi)-sd(diff)/sqrt(5), 4), ";", 
    round(mean(RKHS_CV3$PredAbi-GBLUP_CV3$PredAbi) + (sd(diff)/sqrt(5)), 4),"]") 

# CV4
cat("[", round(mean(RKHS_CV4$PredAbi-GBLUP_CV4$PredAbi)-sd(diff)/sqrt(5), 4), ";", 
    round(mean(RKHS_CV4$PredAbi-GBLUP_CV4$PredAbi) + (sd(diff)/sqrt(5)), 4),"]") 

# CV5
cat("[", round(mean(RKHS_CV5$PredAbi-GBLUP_CV5$PredAbi)-sd(diff)/sqrt(5), 4), ";", 
    round(mean(RKHS_CV5$PredAbi-GBLUP_CV5$PredAbi) + (sd(diff)/sqrt(5)), 4),"]") 

#-------------------------------------------------------------------------------
# 9 - Recovery of entries by rank
#-------------------------------------------------------------------------------
# Summarize the recovery by rank results
recoveryResults <- rbind(
      GBLUP_CV1=rowMeans(colMeans(GBLUP_CV1$topRecovery)),
      GBLUP_CV2=rowMeans(colMeans(GBLUP_CV2$topRecovery)),
      GBLUP_CV3=rowMeans(colMeans(GBLUP_CV3$topRecovery)),
      GBLUP_CV4=rowMeans(colMeans(GBLUP_CV4$topRecovery)),
      GBLUP_CV5=rowMeans(colMeans(GBLUP_CV5$topRecovery)),
      RKHS_CV1=rowMeans(colMeans(RKHS_CV1$topRecovery)),
      RKHS_CV2=rowMeans(colMeans(RKHS_CV2$topRecovery)),
      RKHS_CV3=rowMeans(colMeans(RKHS_CV3$topRecovery)),
      RKHS_CV4=rowMeans(colMeans(RKHS_CV4$topRecovery)),
      RKHS_CV5=rowMeans(colMeans(RKHS_CV5$topRecovery)))

plot(recoveryResults[1,], type="l", col="blue", xaxt="n", ylim=c(0,1), 
     ylab="% Recovered plants", main="% Recovery of entries")
axis(side=1, at=c(1:7), labels=colnames(df))
lines(recoveryResults[2,], type="l", col="red")
lines(recoveryResults[3,], type="l", col="green")
lines(recoveryResults[4,], type="l", col="black")
lines(recoveryResults[5,], type="l", col="orange")
lines(recoveryResults[6,], lty=2, col="blue")
lines(recoveryResults[7,], lty=2, col="red")
lines(recoveryResults[8,], lty=2, col="green")
lines(recoveryResults[9,], lty=2, col="black")
lines(recoveryResults[10,], lty=2, col="orange")

legend("topleft", cex=0.5, legend=c("GBLUP_CV1", "GBLUP_CV2", "GBLUP_CV3", "GBLUP_CV4", 
       "RKHS_CV5", "RKHS_CV1", "RKHS_CV2", "RKHS_CV3","RKHS_CV4", "RKHS_CV5"),
       text.col=c("blue", "red", "green", "black", "orange", "blue", "red", "green", 
               "black", "orange"), lty=c(rep(1,5), rep(2,5)))

#-------------------------------------------------------------------------------
# 10 - Explore how many CV replications are needed
#-------------------------------------------------------------------------------
scheme2 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION",
                         k=5, replication=20, seed=NULL, exclusive=TRUE, 
                         sampling="incompleteTrial",verbose=TRUE)
# Run CV_x replications with x=number of replication
ETA <- NULL
CV_1 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                      G=G,factor="LOCATION", trait="YIELD", nIter=1500, burnIn=250,
                      replications=1)

CV_3 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                      G=G, factor="LOCATION", trait="YIELD", nIter=1500, 
                      burnIn=250, replications=3)

CV_5 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                      G=G, factor="LOCATION", trait="YIELD", nIter=1500, 
                      burnIn=250, replications=5)

CV_10 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                       G=G, factor="LOCATION", trait="YIELD", nIter=1500, 
                       burnIn=250, replications=10)

CV_15 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                       G=G, factor="LOCATION", trait="YIELD", nIter=1500, 
                       burnIn=250, replications=15)

CV_20 <- inferenceBGLR(P, CVscheme=scheme2, modelName="model1", id="GERMPLASM",
                       G=G, factor="LOCATION", trait="YIELD", nIter=1500, 
                       burnIn=250, replications=20)

save(CV_1, file=file.path(outputPath, "CV1.rda"))
save(CV_3, file=file.path(outputPath, "CV3.rda"))
save(CV_5, file=file.path(outputPath, "CV5.rda"))
save(CV_10, file=file.path(outputPath, "CV10.rda"))
save(CV_15, file=file.path(outputPath, "CV15.rda"))
save(CV_20, file=file.path(outputPath, "CV20.rda"))

# Plot the results
B <- 10 # B: number of re-samples
p <- 0.80 # power
sqrt(p*(1-p))/sqrt(B)

load(file.path(outputPath,"CV1.rda"))
load(file.path(outputPath,"CV3.rda"))
load(file.path(outputPath,"CV5.rda"))
load(file.path(outputPath,"CV10.rda"))
load(file.path(outputPath,"CV15.rda"))

means <- apply(CV_3$PredAbi, 2, mean)
sd <- apply(CV_3$PredAbi, 2, function(x) sd(x)/sqrt(length(x)))
mean(CV_1$PredAbi)
sd(CV_3$PredAbi)/sqrt()

boxplot(apply(CV_3$PredAbi, 1, mean, na.rm=TRUE), 
        apply(CV_5$PredAbi, 1, mean, na.rm=TRUE),
        apply(CV_10$PredAbi, 1, mean, na.rm=TRUE),
        apply(CV_15$PredAbi, 1, mean, na.rm=TRUE))

data <- data.frame(replications=c(1,3,5,10,15),
                   mean=c(mean(CV_1$PredAbi, na.rm=TRUE), 
                          mean(CV_3$PredAbi, na.rm=TRUE), 
                          mean(CV_5$PredAbi, na.rm=TRUE), 
                          mean(CV_10$PredAbi, na.rm=TRUE),
                          mean(CV_15$PredAbi, na.rm=TRUE)),
                   SE=c(sd(CV_1$PredAbi, na.rm=TRUE)/sqrt(1), 
                        sd(CV_3$PredAbi, na.rm=TRUE)/sqrt(3), 
                        sd(CV_5$PredAbi, na.rm=TRUE)/sqrt(5), 
                        sd(CV_10$PredAbi, na.rm=TRUE)/sqrt(10),
                        sd(CV_15$PredAbi, na.rm=TRUE)/sqrt(15)))

#===============================================================================
# End of script                                                                #
#===============================================================================