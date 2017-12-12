#===============================================================================
# Example object script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)       #    #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#===============================================================================
#-------------------------------------------------------------------------------
# exampleCV data
#-------------------------------------------------------------------------------
germplasm <- c(paste("PLANT0", rep(1:9), sep=""), paste("PLANT", rep(10:20), sep=""))
data <- data.frame(GERMPLASM= germplasm, LOCATION=rep(1:5, each=20), 
          trait=rnorm(100, 0,1), RANGE=rep(1:4, each=5), ROW=rep(1:5, each=1))
exampleCV <- data[order(data$GERMPLASM, data$LOCATION),]
save(exampleCV, file="exampleCV.rda")
