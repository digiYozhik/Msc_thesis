#===============================================================================
# Analysis script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)             #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#                                                                              #
# Step 2 : - Prepare genotypes                                                 #
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
library(ggplot2)

#-------------------------------------------------------------------------------
# 01 Read genotype data
#-------------------------------------------------------------------------------
data(genotypes)

#-------------------------------------------------------------------------------
# 02 Prepare and clean genotype data
#-------------------------------------------------------------------------------
#******************************************************************************#
# Following steps are performed on the raw genotypes, but not shared because   #
# of confidentiality reasons of functionalities, marker and genetic maps       # 
# information:                                                                 #
# - Conversion of the marker alleles to numeric scores. Conversion was done    # 
#   relative to the minor allele of the markers.                               #
# - Summary of the data                                                        #
#   + Graphical genotypes were checked for abnormalities, but looked ok        #
#   + 570 markers were found with duplicated genotypes (ignoring missingness), # 
#     No entries showed duplication based on genotype information              #
# - Genotype imputation of the markers                                         #
#   + markers were imputed using a hidden Markov Model (HMM) which takes into  #
#     account the information of the cross-type, the order of the markers on   #
#     a genetic map, and the observed recombinations                           #
# - Marker QC                                                                  #
#   + The expected Mendelian segregation ratio's were calculated for a BC3F2   #
#     cross and compared to the  observed ratio's.                             #
#     - Expected segregation ratio: 0:0.90625, 1:0.0625, 2:0.03125             #
#     - Observed segregation ratio: 0:0.90625, 1:0.0625, 2:0.03125             #
#     - A chi-square test was performed using chisq.test:                      #
#       chisq.test(x=c(0.0695,0.084,0.22), y=c(0.90625,0.0625,0.03125)) with   #
#       obtained p-value=0.1991                                                #
#     - Further, an overall of 24% of the markers showed segregation distortion# 
#     - Marker segregation patterns were checked against a genetic map, where  #
#       clearly some regions showed up (which most likely refer to selection   #
#       pressure in the process of crossing). No markers were excluded on the  # 
#       basis of their segregation.                                            #
#   + No markers with more than 10% missing genotype information (NA's) were   #
#      detected. Furhter, an overall missingness of 0.06% was observed.        #
#   + Markers were checked on their minor allele frequencies (shown in below   #
#     code). None were excluded on basis of below set criterium                #
#   + The LD extent was checked but not reported                               #
# - The numerical genotype matrix was exported as Mraw. Users can use          #
#   >data(Mraw) to access the data (">" means console).                        #
#******************************************************************************#

#-------------------------------------------------------------------------------
# 3 - Minor allele frequency (MAF) calculations
#-------------------------------------------------------------------------------
data(Mraw)
# Set a minimum threshold for the MAF
minMAF <- 1/100

# Check if markers need to be remove using minMAF criterium
alleleFrequencies <- rowMeans(Mraw,na.rm=TRUE)/2
MAF <- ifelse(alleleFrequencies > 0.5, 1-alleleFrequencies, alleleFrequencies)
MAFclean <- ifelse(alleleFrequencies > minMAF, 
                   1-alleleFrequencies, alleleFrequencies)
toRemove <- MAFclean < minMAF
sum(toRemove)

# Plot the minor allele frequencies
MAF <- ifelse(alleleFrequencies > 0.5, 1-alleleFrequencies, alleleFrequencies)

png(file.path(outputPath, "MAF.png"))
  MAF <- data.frame(maf=MAF)
  ggplot(MAF, aes(x=maf)) + geom_histogram(colour = "darkgreen", fill = "white", 
    binwidth = 0.01) + scale_x_continuous(lim=c(0,0.5)) + 
  ggtitle('Histogram of Minor Allele Frequency')
dev.off()

#-------------------------------------------------------------------------------
# 4 - Prepare and summarize the genotypic relationship matrix (G-matrix)
#-------------------------------------------------------------------------------
Graw <- makeG(t(Mraw), bAlleleFrequency=rowMeans(Mraw))
dim(Graw) #299 x 299 entries genetic relationships

# Summarize G matrix
## Eigen Vector Decomposition (EVD)
EVD <- eigen(Graw)
varianceExplained <- cumsum(EVD$values)/sum(EVD$values)

plot(EVD$values,main='Eigenvalues',ylab='Eigenvalue',col=2)
plot(varianceExplained,main='Variance Explained',
     ylab='Cumulative Sum',col=2) 
abline(h=c(0.5,0.75,0.9),col=4,lty=2)

png(file.path(outputPath, "EVD_G.png"))
plot(EVD$vectors[,1:2], main='Eigenvectors',
     xlab=paste('1st Eigenvector (%VE:',round(varianceExplained[1],2)*100, "%)"),
     ylab=paste('2nd Eigenvector (%VE:',round(varianceExplained[2]-
          varianceExplained[1],2)*100, "%)"),col=2, pch=19)
dev.off()

## Histograms
hist(Graw)
abline(v=mean(Graw), col="red", lty=2)
mean(Graw)
max(Graw)
min(Graw)

# Frequency distribution for the diagonal element
hist(diag(Graw))
abline(v=mean(diag(Graw)), col="red", lty=2)
mean(diag(Graw))

## Multi-Dimensional Scaling (MDS)
d <- dist(Graw) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- as.numeric(fit$points[,1])
y <- as.numeric(fit$points[,2])

png(file.path(outputPath, "MDS_Gmatrix.png"))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric MDS of G-matrix", type="n")
text(x, y, labels = gsub("LINE", "", row.names(Graw)), cex=0.7)
dev.off()

## Heatmap of the G-matrix
png(file.path(outputPath, "Gmatrix.png"))
heatmap(Graw,Rowv=NA, Colv=NA, labRow=NA, labCol=NA, main="G-matrix",
        revC=TRUE)
dev.off()
save(Graw, file=file.path(dataPath, "Graw.rda"))

#===============================================================================
# End of script                                                                #
#===============================================================================