#' Make a centered (scaled) genetic marker matrix
#' 
#' 
#'@param M matrix of genotypes in bi-allelic format (e.g. AA, AC, CC). Missing 
#'data is entered as NA.
#'@param bAlleleFrequency vector of B allele frequencies for the markers. Default 
#'is NULL, which takes the column means of the M matrix or the column means divided 
#'by two in case of haplotypes.
#'@param scaled logical whether scaling of genotype score need to be done. Default 
#'is TRUE.
#'@param haploid logical wheter marker scores are haplotypes. Default is FALSE.
#'@param SNPinRows logical whether markers are in rows. If FALSE markers are 
#'assumed to be present in columns. Default is TRUE.
#'@return realization of the Z matrix in R matrix format
#'@export
#'@examples
#'data(M)
#'head(makeZ(M))
#'
makeZ <- function(M, bAlleleFrequency=NULL, scaled=TRUE, haploid=FALSE, 
                  SNPinRows=TRUE) {
  if(SNPinRows==TRUE)
    M <- t(M) 
  
  if(is.null(bAlleleFrequency)) { 
    if(haploid==FALSE) { 
      bAlleleFrequency <- colMeans(M)/2 
    } else { 
      bAlleleFrequency <- colMeans(M) 
    }
  }
  Z <- t(apply(M, 1, function(x) {return(x-bAlleleFrequency)}))
  varg <- 2*sum(bAlleleFrequency*(1-bAlleleFrequency))
  if(scaled==TRUE)
    Z <- Z/sqrt(varg)
  
  if(SNPinRows==TRUE) { 
    rownames(Z) <- rownames(M)
    colnames(Z) <- colnames(M)
  } else {
    rownames(Z) <- rownames(M)
    colnames(Z) <- colnames(M)
  }
  return(Z)
}

#' Make genomic relationship matrix using Van Raden (2008)
#' 
#' 
#'@param x matrix of genotypes in bi-allelic format (e.g. AA, AC, CC). Missing 
#'data is entered as NA.
#'@param bAlleleFrequency vector of B allele frequencies for the markers. Default 
#'is NULL which uses the column means of the vector x in the calculation of the 
#'G matrix.
#'@return realization of the G matrix in R matrix format
#'@details The function first calculates the Z matrix using makeZ, then scales 
#'the Z matrix using bAlleleFrequency to obtain the G matrix as defined by 
#'Van Raden (2008).
#'@export
#'@seealso makeZ
#'@examples
#'data(M)
#'head(makeG(M))
#'@references 
#'Van Raden, P.M. (2008). Efficient methods to compute genomic predictions. 
#'Journal of Dairy Science, 91(11): 4412-4423.
#' 
makeG <- function (x, bAlleleFrequency=NULL) {
  nSNP <- ncol(x) 
  if (!is.null(bAlleleFrequency) && length(bAlleleFrequency) != nSNP) 
    stop("maf vector does not have correct dimensions")
  #  if (SNPinRows==TRUE) {  x <- t(x) }
  if (is.null(bAlleleFrequency)) {
    bAlleleFrequency <- colMeans(x)/2 }
  Z <- t(apply(x, 1, function(y) {
    return(y - 2 * bAlleleFrequency)
  }))
  G <- tcrossprod(Z)/(2 * sum(bAlleleFrequency * (1 - bAlleleFrequency)))
  return(G)
}
