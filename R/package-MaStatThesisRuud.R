#' Introductory comments on MaStatThesisRuud
#'
#' @details The package holds the data and user defined functionalities used in 
#' the analyses in the M.Sc. thesis of Ruud Derijcker (see reference 3). 
#' 
#' Use the \code{example} function to run examples of the various functions. 
#' To view the code of a function simply type the function name in the console.
#' @keywords GS BGLR
#' @name MaStatThesisRuud
#' @docType package
#' @title MaStatThesisRuud
#' @import synbreed
#' @import coda
#' @import BGLR
#' @import ggplot2
#' @author Ruud Derijcker <ruud.derijcker@@ugent.be>
#' @references \describe{
#'   \item{\code{1}:}{De Los Campos, G., Perez, P.(2014). BGLR: Bayesian 
#'   Generalized Linear Regression. Version 1.0.3. 
#'   (http://CRAN.R-project.org/package=BGLR).}
#'   \item{\code{2}:} {Jarquin, D. et al. (2014). A reaction norm model for 
#'   genomic selection using high-dimensional genomic and environmental data. 
#'   Theoretical and Applied Genetics 127(3):595-607.}
#'   \item{\code{3}:} {Derijcker, R. (2015). Investigating incorporation of 
#'   genotype x environment interaction (G x E) for genomic selection in a 
#'   practical setting. Unpublished M.Sc. thesis. University of Ghent: Belgium.}    
#'}
NULL

#' Raw phenotype data used in thesis
#' 
#' Raw phenotype data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi--parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Yield performance 
#' measurements were obtained per plot. The data frame contains the field design 
#' information and the yield performance values averaged per plot. The data frame  
#' holds information of 2310 observations and 9 features. The features are 
#' detailed below and represent the columns in the data frame.
#' \describe{
#'      \item{\code{GERMPLASM}:}{The entry names.}
#'      \item{\code{LOCATION}:}{The name of the locations.}
#'      \item{\code{RANGE}:}{The range coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{ROW}:}{The row coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{RANGEROW}:}{A combination of the range and row coordinates.}
#'      \item{\code{LOCAL_ROW}:}{The coordinates for the rows linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{LOCAL_RANGE}:}{The coordinates for the ranges linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{PLOT}:}{The reference to the plot of the observation.}
#'      \item{\code{YIELD}:}{The average yield performance measures of the plots for 
#'      the respective observations.}
#' }
#' @author Ruud Derijcker
#' @name data-phenotypes
#' @docType data
#' @usage phenotypes
#' @keywords phenotypes
#' @examples 
#' data(phenotypes)
#' head(phenotypes)
#' 
NULL

#' Modeled phenotypic data used in thesis
#' 
#' Modeled phenotypic data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi--parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Yield performance 
#' measurements were averaged per plot. Location 7 was excluded from the 
#' analysis. Next, a linear mixed model was applied to obtain BLUP predictions 
#' for the yield performance measures where the yield performance measures were 
#' adjusted using the information of the checks, and where a spatial (AR1xAR1) 
#' covariance structure was applied for the rows and ranges in the fields. 
#' The resulting data frame holds information of 1774 observations and 19 features.
#' The features are detailed below and represent the columns in the data frame.
#' \describe{
#'      \item{\code{MERGE}:}{Observation names, which are the combined names of 
#'      the entries and the locations.}
#'      \item{\code{GERMPLASM}:}{The entry names.}
#'      \item{\code{LOCATION}:}{The name of the locations.}
#'      \item{\code{OVERALL_RANGE}:}{The range coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{OVERALL_ROW}:}{The row coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{RANGEROW}:}{A combination of the range and row coordinates.}
#'      \item{\code{ROW}:}{The coordinates for the rows linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{RANGE}:}{The coordinates for the ranges linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{PLOT}:}{The reference to the plot of the observation.}
#'      \item{\code{YIELD}:}{The average yield performance measures of the plots for 
#'      the respective observations.}
#'      \item{\code{EXPERIMENT}:}{Name of the experiment.}
#'      \item{\code{SUBSUBBLOCK}:}{Placeholder in case subsubblocks were defined 
#'      in case of blocking. Not used here.}
#'      \item{\code{SUBBLOCK}:}{Placeholder in case subblocks were defined 
#'      in case of blocking. Not used here.}
#'      \item{\code{BLOCK}:}{Placeholder in case blocks were defined in case of 
#'      blocking. Not used here.}
#'      \item{\code{CHECK}:}{Numeric showing whether the observation is a check. 
#'      Here the checks were already filtered out.}
#'      \item{\code{ROWf}:}{Factor for ROW.}
#'      \item{\code{RANGEf}:}{Factor for RANGE.}
#'      \item{\code{residuals}:}{Residuals obtained after applying the above 
#'      described linear mixed model.}
#'      \item{\code{predicted.values}:}{BLUP predictions from the fitted above 
#'      described linear mixed model.}
#' }
#' @author Ruud Derijcker
#' @name data-Praw
#' @docType data
#' @usage Praw
#' @keywords Praw
#' @examples 
#' data(Praw)
#' str(Praw)
#' 
NULL

#' Modeled cleaned phenotypic data used in thesis
#' 
#' Modeled cleaned phenotypic data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi--parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Yield performance 
#' measurements were averaged per plot. Location 7 was excluded from the 
#' analysis. Next, a linear mixed model was applied to obtain BLUP predictions 
#' for the yield performance measures where the yield performance measures were 
#' adjusted using the information of the checks, and where a spatial (AR1xAR1) 
#' covariance structure was applied for the rows and ranges in the fields. After 
#' aligning the phenotypes with the genotypes, the resulting data frame holds 
#' information of 1768 observations and 19 features. The features are detailed 
#' below and represent the columns in the data frame.
#' \describe{
#'      \item{\code{MERGE}:}{Observation names, which are the combined names of 
#'      the entries and the locations.}
#'      \item{\code{GERMPLASM}:}{The entry names.}
#'      \item{\code{LOCATION}:}{The name of the locations.}
#'      \item{\code{OVERALL_RANGE}:}{The range coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{OVERALL_ROW}:}{The row coordinates when all fields are 
#'      seen as being part of one big field, i.e. same reference grid for all fields.}
#'      \item{\code{RANGEROW}:}{A combination of the range and row coordinates.}
#'      \item{\code{ROW}:}{The coordinates for the rows linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{RANGE}:}{The coordinates for the ranges linked to the locations. 
#'      Here the reference grid is the location itself.}
#'      \item{\code{PLOT}:}{The reference to the plot of the observation.}
#'      \item{\code{YIELD}:}{The average yield performance measures of the plots for 
#'      the respective observations.}
#'      \item{\code{EXPERIMENT}:}{Name of the experiment.}
#'      \item{\code{SUBSUBBLOCK}:}{Placeholder in case subsubblocks were defined 
#'      in case of blocking. Not used here.}
#'      \item{\code{SUBBLOCK}:}{Placeholder in case subblocks were defined 
#'      in case of blocking. Not used here.}
#'      \item{\code{BLOCK}:}{Placeholder in case blocks were defined in case of 
#'      blocking. Not used here.}
#'      \item{\code{CHECK}:}{Numeric showing whether the observation is a check. 
#'      Here the checks were already filtered out.}
#'      \item{\code{ROWf}:}{Factor for ROW.}
#'      \item{\code{RANGEf}:}{Factor for RANGE.}
#'      \item{\code{residuals}:}{Residuals obtained after applying the above 
#'      described linear mixed model.}
#'      \item{\code{predicted.values}:}{BLUP predictions from the fitted above 
#'      described linear mixed model.}
#' }
#' @author Ruud Derijcker
#' @name data-P
#' @docType data
#' @usage P
#' @keywords P
#' @examples 
#' data(P)
#' str(P)
#' 
NULL

#' Raw SNP marker data used in thesis
#' 
#' Raw SNP marker data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi-parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Genotypes were obtained 
#' for the entries on one locations from a set of genome-wide spread SNP 
#' markers which were scored in a bi-allelic way and represented as a combination 
#' of two of the four bases that make up the DNA (e.g. AA, AC, CC). The data frame 
#' holds information about 1088 SNP markers (rows) for 305 entries (columns), 
#' including checks entries.
#' @author Ruud Derijcker
#' @name data-genotypes
#' @docType data
#' @usage genotypes
#' @keywords genotypes
#' @examples 
#' data(genotypes)
#' head(genotypes)
#' 
NULL

#' Numeric SNP marker data used in thesis
#' 
#' Numeric SNP marker data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi-parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Genotypes were obtained 
#' for the entries on one locations from a set of genome--wide spread SNP 
#' markers which were scored in a bi-allelic way. After quality control, 
#' genotypes for 299 cross-entries were obtained for 1088 SNP markers. 
#' The observed genotype matrix is converted to allelic content by counting 
#' the number of B alleles for every cell in the genotype matrix. The obtained 
#' matrix of {0,1,2}-counts is called the M matrix.
#' @author Ruud Derijcker
#' @name data-Mraw
#' @docType data
#' @usage Mraw
#' @keywords Mraw
#' @examples 
#' data(Mraw)
#' str(Mraw)
#' 
NULL

#' Cleaned numeric SNP marker data used in thesis
#' 
#' Cleaned numeric SNP marker data used in thesis
#' 
#' @details The data comes from a field trial that is part of a cotton 
#' breeding program. The trial was set up in 2012 across 7 locations in the 
#' US Cotton Belt. At every location the same bi-parental BC_3F_2 was grown 
#' together with a number of entries serving as checks. Genotypes were obtained 
#' for the entries on one locations from a set of genome--wide spread SNP 
#' markers which were scored in a bi-allelic way. After quality control, 
#' genotypes for 299 cross-entries were obtained for 1088 SNP markers. 
#' The observed genotype matrix is converted to allelic content by counting 
#' the number of B alleles for every cell in the genotype matrix. When aligning 
#' the phenotype and genotype information with respect to the entries, we reduce 
#' the number of entries to 296 by discarding genotyped entries that have 
#' completely unobserved phenotypes.
#' @author Ruud Derijcker
#' @name data-M
#' @docType data
#' @usage M
#' @keywords M
#' @examples 
#' data(M)
#' str(M)
#' 
NULL

#' Realization of the G-matrix used in thesis
#' 
#' Realization of the G-matrix used in thesis
#' 
#' @details The realization of the G-matrix was calculated using the numeric 
#' SNP genotype matrix M following the calculation detailed by Van Raden (2008). 
#' The dimension are 296 x 296 where both rows and columns represent the entries. 
#' The elements of the matrix show the genetic relationships between the entries 
#' that are a combination of a row and column coordinate. 
#' @author Ruud Derijcker
#' @name data-G
#' @docType data
#' @usage G
#' @keywords G
#' @seealso makeG
#' @examples 
#' data(G)
#' str(G)
#' 
NULL

#' Example data for CV functionalities used in thesis
#' 
#' Example data for CV functionalities used in thesis
#' 
#' @details The example data frame contains 100 entries grown at 5 locations such 
#' that 20 entries are grown per location. The trait was obtained using the 
#' rnorm function. The data frame holds information for the 100 entries, and 5 
#' features represented in the columns and detailed below:
#' \describe{
#'      \item{\code{GERMPLASM}:}{The 5 x 20 entry names.}
#'      \item{\code{LOCATION}:}{Numeric describing which of the 5 locations the 
#'      entry is coming fromm.}
#'      \item{\code{trait}:}{A value for the trait obtained by rnorm(100, 0, 1).}
#'      \item{\code{RANGE}:}{The coordinates for the ranges linked to the locations.}
#'      \item{\code{ROW}:}{The coordinates for the rows linked to the locations.}
#' }
#' @author Ruud Derijcker
#' @name data-exampleCV
#' @docType data
#' @usage exampleCV
#' @keywords exampleCV
#' @examples 
#' data(exampleCV)
#' head(exampleCV)
#' 
NULL
