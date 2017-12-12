#===============================================================================
# Installation script for M.Sc. thesis of Ruud Derijcker (2014 - 2015)         #    #
# "Investigating incorporation of genotype x environment interaction (G x E)   #
#  for genomic selection in a practical setting."                              #
#===============================================================================
library(devtools)
pkgPath <- "/home/bbrud/ThesisPackage/pkg"
downloadPath <- "/home/bbrud/ThesisPackage/etc/downloads"

# Check package structure
check(pkg=pkgPath, check_dir=downloadPath, cleanup=FALSE, cran=FALSE)

# Build Package Source
build(pkg=pkgPath, binary=TRUE, vignettes=FALSE, path=downloadPath)

# Build Package Source (no binary)
build(pkg=pkgPath, binary=FALSE, vignettes=FALSE, path=downloadPath)
