#===============================================================================
# Analysis script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)             #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#                                                                              #
# Step 1 : - Exploration of experimental design                                #
#          - Processing of phenotype information                               #
#          - Linear mixed modeling (spatial correction, BLUP predictions,      #
#            heritability calculations)                                        #
#          - G x E exploration                                                 #
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
# David Butler (2009). myf: Utility functions for asreml objects.
# R package version 1.0.
library(myf) # for G x E inspection 
# library(asreml)

#-------------------------------------------------------------------------------
# 1 - Import phenotype data
#-------------------------------------------------------------------------------
data(phenotypes) # raw phenotype data stored in MaStatThesis package
# Format the raw data
rawData <- phenotypes
colnames(rawData)[match("RANGE", colnames(rawData))] <- "OVERALL_RANGE"
colnames(rawData)[match("ROW", colnames(rawData))] <- "OVERALL_ROW"
colnames(rawData)[match("LOCAL_RANGE", colnames(rawData))] <- "RANGE"
colnames(rawData)[match("LOCAL_ROW", colnames(rawData))] <- "ROW"
rawData$EXPT <- "THESIS"
rawData$BLOCK <- rawData$SUBBLOCK <- rawData$SUBSUBBLOCK <- NA
rawData$LOCATION <- as.character(rawData$LOCATION)
checks <- c("CHECK01","CHECK02","CHECK03","CHECK04")
rawData$CHECK <- ifelse(rawData$GERMPLASM=="CHECK01", 1, 
                        ifelse(rawData$GERMPLASM=="CHECK02", 2, 
                        ifelse(rawData$GERMPLASM=="CHECK03", 3,
                        ifelse(rawData$GERMPLASM=="CHECK04", 4,0))))

#-------------------------------------------------------------------------------
# 2 - Data exploration
#-------------------------------------------------------------------------------
# Data consistency checks
table(rawData$GERMPLASM, rawData$LOCATION)
table(rawData$LOCATION)
table(rawData$GERMPLASM)
table(rawData$CHECK)
summary(rawData$YIELD)
fields <- as.character(unique(rawData$LOCATION)); fields
table(rawData$CHECK)

# Exploration
## Fields design plots
pdf(file=file.path(outputPath, "Fieldplots.pdf"), 15, 8)
for(iField in 1: length(fields)){
  temp <- rawData[rawData$LOCATION %in% fields[iField],]
  temp$CHECK <- as.factor(temp$CHECK)
  temp$GERMPLASM <- gsub("LINE", "E", temp$GERMPLASM)
  temp$GERMPLASM <- gsub("CHECK", "C", temp$GERMPLASM)
  print(plotFields(x=temp, column="CHECK", textColumn="GERMPLASM",
             main=paste("Fieldplot for location", fields[iField]),
             colors=c("grey", "red","blue","green","yellow")))
}
dev.off()

## Boxplots for yield performance including checks alongside boxes
png(file.path(outputPath, "boxplots_allChecks.png"), width = 1200, height = 800)
trait <- "YIELD"
boxplots_with_checks(x=rawData, traits="YIELD", factor = "LOCATION", 
                     checks=rawData$GERMPLASM[rawData$CHECK!=0], 
                     checks_factor = "GERMPLASM", checks.color = "grey", 
                     boxcolor="blue", main="Yield by location\n(all checks)")
dev.off()

# boxplots for yield with only CHECK04
png(file.path(outputPath, "boxplots_check4.png"), width = 1200, height = 800)
trait <- "YIELD"
boxplots_with_checks(x=rawData, trait="YIELD", factor = "LOCATION", 
                     checks=rawData$GERMPLASM[rawData$CHECK==4], 
                     checks_factor = "GERMPLASM", checks.color = "grey", 
                     boxcolor = "blue", main="Yield by location\n(CHECK04)")
dev.off()

# check unobserved yield data
unobservedYieldData <- vector()
for(iField in 1:length(fields)){
  data <- subset(rawData, LOCATION==fields[iField])
  unobservedYieldData[iField] <- round(sum(is.na(data$YIELD))/nrow(data)*100,2)
}
unobservedYieldData # LOCATION 7 has 73% missing data, leave out in analysis

#-------------------------------------------------------------------------------
# 3 - Subset the data
#-------------------------------------------------------------------------------
rawData <- rawData[rawData$LOCATION != "LOCATION7", ] # 1980 x 14

#-------------------------------------------------------------------------------
# 4 - Check if we benefit from applying a spatial model to the data
#-------------------------------------------------------------------------------
# Prepare data
rawData$ROWf <- as.factor(rawData$ROW)
rawData$RANGEf <- as.factor(rawData$RANGE)
rawData <- rawData[order(rawData$LOCATION, rawData$ROWf, rawData$RANGEf), ]
rawData <- droplevels(rawData)

# Fit a linear mixed model using ASReml-R
require(asreml)
model <- asreml(fixed = YIELD ~ 1,
                random = ~ diag(LOCATION):GERMPLASM,  
                rcov= ~ at(LOCATION):ar1(ROWf):ar1(RANGEf), 
                data= rawData,
                control = asreml.control(workspace=40e6, maxiter=500),
                na.method.X = "include", na.method.Y="include")
# LRT test
# The likelihood ratio test is performed by comparing two times the difference
# in log-likelihood with the quantiles if the \eqn{\chi^2}-distribution with
# as degrees of freedom the difference in number of variance components
# estimated in both models.

## Model with no spatial correction
model1 <- asreml(fixed = YIELD ~ 1,
                random = ~ diag(LOCATION):GERMPLASM,  
                data= rawData,
                control = asreml.control(workspace=40e6, maxiter=500),
                na.method.X = "include", na.method.Y="include")

# Model with spatial correction
model2 <- asreml(fixed = YIELD ~ 1,
                random = ~ diag(LOCATION):GERMPLASM,  
                rcov= ~ at(LOCATION):ar1(ROWf):ar1(RANGEf), 
                data= rawData,
                control = asreml.control(workspace=40e6, maxiter=500),
                na.method.X = "include", na.method.Y="include")
LR.test.asreml(model1, model2)
# > LR.test.asreml(model1, model2)
# 
# Likelihood ratio test (Chi-square)
# 
# data:  model1 and model2
# X = 926.41, degrees of freedom = 17, p-value < 2.2e-16
# alternative hypothesis: Extended model outperforms the 'null model'.

# Conclusion: We go for a spatial model!

#-------------------------------------------------------------------------------
# 5 - Apply a linear mixed model on the data to predict phenotypes and calculate 
#     BLUPs for the yield observations
#-------------------------------------------------------------------------------
# We want to adjust using:
# - check information (i.e. replication)
# - assumed spatial AR1(rows) x AR1(ranges) trend in the field
require(asreml)
model <- asreml(fixed = YIELD ~ 1,
                random = ~ diag(LOCATION):GERMPLASM,  
                rcov= ~ at(LOCATION):ar1(ROWf):ar1(RANGEf), 
                data= rawData,
                control = asreml.control(workspace=40e6, maxiter=500),
                na.method.X = "include", na.method.Y="include")

adjustedYield <- predict(model, classify="GERMPLASM:LOCATION")$predictions$pvals[,1:4]
table(adjustedYield$GERMPLASM, adjustedYield$LOCATION)

# compare BLUPs to crude averages
crudeAverages <- aggregate( x=rawData$YIELD, by= list(GERMPLASM = rawData$GERMPLASM, 
                    LOCATION = rawData$LOCATION), FUN=mean, na.rm=TRUE)
colnames(crudeAverages)[3] <- "crudeAverage"
Praw <- merge(crudeAverages, adjustedYield)
Praw <- Praw[!is.na(Praw$crudeAverage),]
Praw
# > head(Praw)
# GERMPLASM  LOCATION crudeAverage predicted.value standard.error
# 1   CHECK01 LOCATION1      1127.67          1127.1         74.873
# 2   CHECK01 LOCATION2      1806.67          1565.2         60.883
# 3   CHECK01 LOCATION3      1535.83          1476.1         39.866
# 4   CHECK01 LOCATION4       951.33          1354.6         30.192
# 5   CHECK01 LOCATION5      1268.00          1451.9         23.347
# 6   CHECK01 LOCATION6      1422.67          1333.0         83.541

#-------------------------------------------------------------------------------
# 6 - Summarize based on BLUPs
#-------------------------------------------------------------------------------
#******************************************************************************#
# Users can use > data(Praw) in case no asreml license (">" means console)     #
#******************************************************************************#
# Plot the correlation between crude and BLUP averages across locations
png(file.path(outputPath, "adjustedPhenotypes.png"))
plot(Praw$crudeAverage, Praw$predicted.value, xlab="Raw phenotype values (lb/ac)",
     ylab="Predicted phenotype values (lb/ac)", 
     main=paste("Correlation between crude and BLUP averages across locations",
                "\ncorrelation=",round(cor(Praw$crudeAverage, Praw$predicted.value),2),
                sep=""))
dev.off()
# Found correlation: 0.68

# Plot histogram and qqplot of the adjusted yield values
png(file.path(outputPath, "histogramYield.png" ))
par(mfrow=c(1,2))
  hist(Praw$predicted.value, las=2, main="Standardized yield", 
       xlab="Adjusted yield", ylab="Frequency")
  check04 <- mean(Praw$predicted.value[which(Praw$GERMPLASM=="CHECK04")])
  checkOverall <- mean(Praw$predicted.value[which(Praw$GERMPLASM %in% checks)])
  abline(v=check04, col="red")
  abline(v=checkOverall, col="blue")
  legend("topleft", legend=c("CHECK04","All checks"), text.col=c("red","blue"), 
         bty="n")
  qqnorm(Praw$predicted.value, main="Q-Q plot for adjusted yield")
  qqline(Praw$predicted.value)
dev.off()

# Plot the line abundance
png(file.path(outputPath, "lineAbundance.png" ))
plot(table(Praw$LOCATION),xlab="",las=2,ylab='Number of records',
     main="Number of plants per location")
dev.off()

# Plot the correlation in yield between locations of the predicted yield values
PWide <- reshape(Praw[, c("GERMPLASM", "LOCATION", "predicted.value")], 
                 timevar="LOCATION", idvar="GERMPLASM", direction="wide")
png(file.path(texFigurePath, "corplot.png"))
  corplot.data.frame(PWide[, -1], main ="Correlation of yield (lb/ac) by location")
dev.off()

# Store the predicted values
Praw$MERGE <- paste(Praw$GERMPLASM, Praw$LOCATION, sep="_")
Praw <- Praw[,c("predicted.value","MERGE")]
rawData$MERGE <- paste(rawData$GERMPLASM, rawData$LOCATION, sep="_")
Praw <- merge(rawData, Praw, by="MERGE")
Praw <- Praw[!(Praw$GERMPLASM %in% checks),]
save(Praw, file=file.path(outputPath, "Praw.rda"))

# Check the precision of the phenotypes
## Before spatial adjustments (raw data)
FACTOR <- as.factor(as.character(rawData[, "LOCATION"]))
temp <- data.frame(Factor = levels(FACTOR))
means <- std <- cv <- temp
means <- data.frame(means, aggregate(rawData[, "YIELD"], 
                                     list(FACTOR), FUN = mean, na.rm = TRUE)$x)
colnames(means)[ncol(means)] <- "YIELD"
std <- data.frame(std, aggregate(rawData[, "YIELD"], 
                                 list(FACTOR), FUN = sd, na.rm = TRUE)$x)
colnames(std)[ncol(std)] <- "YIELD"
cv <- data.frame(as.character(std$Factor), 100*std$YIELD/means$YIELD)
colnames(cv)[ncol(cv)] <- "YIELD"
colnames(cv)[1] <- "Factor"
rownames(means) <- paste("Mean", rownames(means), sep = ".")
rownames(std) <- paste("SD", rownames(std), sep = ".")
rownames(cv) <-  paste("cv", rownames(cv), sep = ".")
resultBeforeModel <- cbind(means, std$YIELD, cv$YIELD)
rownames(resultBeforeModel) <- NULL
resultBeforeModel

## After spatial adjustment
FACTOR <- as.factor(as.character(Praw[, "LOCATION"]))
temp <- data.frame(Factor = levels(FACTOR))
means <- std <- cv <- temp
means <- data.frame(means, aggregate(Praw[, "predicted.value"], 
                                     list(FACTOR), FUN = mean, na.rm = TRUE)$x)
colnames(means)[ncol(means)] <- "predicted.value"
std <- data.frame(std, aggregate(Praw[, "predicted.value"], 
                                 list(FACTOR), FUN = sd, na.rm = TRUE)$x)
colnames(std)[ncol(std)] <- "predicted.value"
cv <- data.frame(as.character(std$Factor), 100*std$predicted.value/means$predicted.value)
colnames(cv)[ncol(cv)] <- "predicted.value"
colnames(cv)[1] <- "Factor"
rownames(means) <- paste("Mean", rownames(means), sep = ".")
rownames(std) <- paste("SD", rownames(std), sep = ".")
rownames(cv) <-  paste("cv", rownames(cv), sep = ".")
resultAfterModel <- cbind(means, std$predicted.value, cv$predicted.value)
rownames(resultAfterModel) <- NULL
resultAfterModel

## Plot %cv
png(file.path(outputPath, "cv.png"))
  boxcolors <- rainbow(length(fields))
  par(mfcol = c(1, 1), mar = c(8.5, 4, 4, 2))
  v <- summary(model)$varcomp
  v$name <- rownames(v)
  sd <- sqrt(v[grep("!variance", v$name), 2])
  means <- resultAfterModel$predicted.value
  cv <- 100*sd/means
  b <- barplot(cv, col = boxcolors, las = 2, ylim = c(0, 50), 
               main = "Residual coefficient of variation (%)")
  axis(side=1, at=b, labels=levels(as.factor(fields[-length(fields)])), las=2)
dev.off()

#-------------------------------------------------------------------------------
# 7 - Model diagnostics
#-------------------------------------------------------------------------------
rawData$residuals <- model$residuals
help(plotResiduals)
# Make diagnostic plots
pdf(file.path(outputPath, "residualPlots.pdf"))
par(mfcol = c(3, 2))
plotResiduals(x=rawData, residuals="residuals", by="LOCATION")
dev.off()

#-------------------------------------------------------------------------------
# 8 - Calculation of heritability
#-------------------------------------------------------------------------------
rawData$GROUP <- rawData$GERMPLASM
n.checks <- length(checks)
for(i in 1:n.checks) 
  rawData$CHECK[rawData$GERMPLASM == checks[i]] <- i
rawData$CHECK[!rawData$GERMPLASM %in% checks] <- n.checks + 1
rawData$CHECK <- factor(rawData$CHECK)

rawData$fRange <- factor(rawData$RANGE)
rawData$fRow <- factor(rawData$ROW)
rawData <- rawData[order(rawData$LOCATION, rawData$fRange, rawData$fRow),]

# Explicit correction for spatial effects
modelSpatial <- asreml(fixed = YIELD ~ 1,
                       random = ~ LOCATION + GERMPLASM,
                       rcov   = ~ at(LOCATION):ar1(fRange):ar1(fRow),
                       data   =  rawData, 
                       na.method.X = "include")
calculateH2(modelSpatial, data = rawData, checks=checks)
# Method STD.Gen STD.GxE R.GxE STD.Res R.Res   H2
# 1   overall  132.01      NA    NA  206.36    NA 0.29
# 2 line-mean  132.01      NA    NA  206.36     6 0.71

#-------------------------------------------------------------------------------
# 9 - Investigate G x E using factor analytic model
#-------------------------------------------------------------------------------
require(myf)
rawData$LOCATION <- as.factor(rawData$LOCATION)
rawData <- rawData[order(rawData$LOCATION, rawData$ROWf, rawData$RANGEf),]

# Fit a factor analytic model using asReml. We only show the model for 4 factors.
METmodelFA4 <- asreml(fixed = YIELD ~ 1,
                       random = ~ fa(LOCATION, 4):GERMPLASM,
                       rcov   = ~ at(LOCATION):ar1(ROWf):ar1(RANGEf),
                       data   =  rawData, 
                       na.method.X = "include", na.method.Y = "include",
                       control = asreml.control(workspace=40e6, maxiter=500))
METmodelFA4 <- update(METmodelFA4) # update needed for convergence

# Summarize FA4 model using myf's summary.fa function
pdf(file.path(outputPath, "fa4Model.pdf"))
summaryMETmodelFA4 <- print(summary.fa(METmodelFA4, trunc.char=NULL, 
                      heatmap.ord='cluster', uniplot=FALSE))
dev.off()

# Print the covariance matrix for
png("faModel.png", width=800, height=800)
summaryMETmodelFA4 <- print(summary.fa(METmodelFA4, trunc.char=NULL, 
                                       heatmap.ord='cluster', uniplot=FALSE))
dev.off()

# Summarize inference
## genetic variances
diag(summaryMETmodelFA4$gammas$fa$Gmat)
# LOCATION1 LOCATION2 LOCATION3 LOCATION4 LOCATION5 LOCATION6 
# 44714     71605     19620     26440     25140     35854 

## Total variance explained by multiplicative part of model
cat("Total variance explained by the model:", 
    round(summaryMETmodelFA4$gammas[[1]]$'total %vaf',2),"%", "\n")
# Total variance explained by the model: 81.43 % 

## Site specific variance explained
cat("Site specific variance explained:", "\n")
summaryMETmodelFA4$gammas[[1]]$'site %vaf'
# fac_1      fac_2     fac_3    fac_4     all
# LOCATION1 26.439 6.5986e+00  0.012756  0.14741  33.198
# LOCATION2 87.838 9.3353e+00  2.327919  0.49853 100.000
# LOCATION3 38.702 5.5637e-05  1.439002  0.85272  40.994
# LOCATION4 64.330 1.8318e+00 17.050684 16.78771 100.000
# LOCATION5 61.209 5.4951e+00 21.952740 11.34348 100.000
# LOCATION6 71.302 1.8012e+01  9.757385  0.92913 100.000

#===============================================================================
# End of script                                                                #
#===============================================================================