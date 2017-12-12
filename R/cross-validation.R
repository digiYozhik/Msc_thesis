#' Function to make a CV scheme based on random sampling of observation IDs
#' 
#' Function to make a CV scheme based on random sampling of observation IDs
#' 
#'@param ID character vector of the observation IDs used in the randomization. 
#'The names are a combination of the entry name and the location.
#'@param seed numeric value for the seed value used for the randomization by the 
#'set.seed function. In this way randomization can be reproduced by the user.
#'Default is NULL, which uses 123 as value for the seed.
#'@param k integer value for the number of folds used in the k-cross-validation.
#'@param exclusive logical whether sampling should be done with replacement. The 
#'argument is passed to the replace argument of the samp.int function as the negation, 
#'i.e. exclusive is TRUE means replace=FALSE, such that the probability of 
#'choosing the next item is proportional to the weights amongst the remaining 
#'items. 
#'@return named vector of numeric scores showing the assignment of the observations 
#'to their respective set used in the k-fold cross-validation.
#'@export
#'@details for the randomization the sample function is used.
#'@author Ruud Derijcker
#'@references Based on synbreed's crossVal function
#'@examples
#'data(exampleCV)
#'y <- exampleCV[,which(colnames(exampleCV) %in% c("GERMPLASM", "LOCATION"))]
#'colnames(y) <- c("IDUnique","FACTOR")
#'y$ID <- paste(y$IDUnique, y$FACTOR, sep="_")
#'y <- na.omit(y)
#'n <- length(y$ID)
#'output <- CVrandomByID(y$ID, seed=123, k=5, exclusive=TRUE)
#'table(output)
#'head(output)
#'
CVrandomByID <- function(ID, seed, k, exclusive) {
  set.seed(seed)
  n <- length(ID)
  modu <- n %% k
  scheme <- sample(c(rep(1:k, each = (n - modu)/k), 
                     sample(1:k, modu)), n, replace = !exclusive)
  names(scheme) <- ID
  return(scheme)
}

#' Function to make a CV scheme based on random sampling of entry IDs 
#' across a specified factor.
#' 
#' Function to make a CV scheme based on random sampling of entry IDs 
#' across a specified factor.
#'  
#'@param ID character vector of the observation IDs. 
#'The names are a combination of the entry name and the location.
#'@param factorID character vector of the entry IDs used in the randomization.
#'@param seed numeric value for the seed value used for the randomization by the 
#'set.seed function. In this way randomization can be reproduced by the user.
#'Default is NULL, which uses 123 as value for the seed.
#'@param k integer value for the number of folds used in the k-cross-validation.
#'@param exclusive logical whether sampling should be done with replacement. The 
#'argument is passed to the replace argument of the samp.int function as the negation, 
#'i.e. exclusive is TRUE means replace=FALSE, such that the probability of 
#'choosing the next item is proportional to the weights amongst the remaining 
#'items. 
#'@return named vector of numeric scores showing the assignment of the observations 
#'to their respective set used in the k-fold cross-validation.
#'@export
#'@details for the randomization the sample function is used.
#'@author Ruud Derijcker
#'@references Based on synbreed's crossVal function
#'@examples
#'data(exampleCV)
#'y <- exampleCV[,which(colnames(exampleCV) %in% c("GERMPLASM", "LOCATION"))]
#'colnames(y) <- c("IDUnique","FACTOR")
#'y$ID <- paste(y$IDUnique, y$FACTOR, sep="_")
#'y <- na.omit(y)
#'n <- length(y$ID)
#'output <- CVrandomAccrossFactor(y$ID, factorID=y$IDUnique, seed=123, 
#'k=5, exclusive=TRUE)
#'table(output)
#'head(output)
#'
CVrandomAccrossFactor <- function(ID, factorID, seed, k, exclusive) {
  set.seed(seed)
  which.pop <- unique(factorID)
  if (length(which.pop) < k) 
    stop("The parameter k can not be greater than the number of populations!")
  y.u <- unique(ID)
  y2 <- matrix(y.u[order(factorID)], ncol = 1)
  b <- table(factorID)
  modu <- length(which.pop)%%k
  val.samp <- sample(c(rep(1:k, each = (length(which.pop) - 
                                          modu)/k), sample(1:k, modu)), length(which.pop), 
                     replace = !exclusive)
  val.samp2 <- rep(val.samp, b)
  val.samp3 <- data.frame(y2, val.samp2)
  scheme <- val.samp3[order(as.character(val.samp3[, 1])), "val.samp2"]
  names(scheme) <- ID
  return(scheme)
}

#' Function to make a CV scheme based on random sampling of entry IDs 
#' within a specified factor.
#'  
#' Function to make a CV scheme based on random sampling of entry IDs 
#' within a specified factor.
#' 
#'@param ID character vector of the observation IDs. 
#'The names are a combination of the entry name and the location.
#'@param factorID character vector of the entry IDs used in the randomization.
#'@param seed numeric value for the seed value used for the randomization by the 
#'set.seed function. In this way randomization can be reproduced by the user.
#'Default is NULL, which uses 123 as value for the seed.
#'@param k integer value for the number of folds used in the k-cross-validation.
#'@param exclusive logical whether sampling should be done with replacement. The 
#'argument is passed to the replace argument of the samp.int function as the negation, 
#'i.e. exclusive is TRUE means replace=FALSE, such that the probability of 
#'choosing the next item is proportional to the weights amongst the remaining 
#'items. 
#'@return named vector of numeric scores showing the assignment of the observations 
#'to their respective set used in the k-fold cross-validation.
#'@export
#'@details for the randomization the sample function is used.
#'@author Ruud Derijcker
#'@references Based on synbreed's crossVal function
#'@examples
#'data(exampleCV)
#'y <- exampleCV[,which(colnames(exampleCV) %in% c("GERMPLASM", "LOCATION"))]
#'colnames(y) <- c("IDUnique","FACTOR")
#'y$ID <- paste(y$IDUnique, y$FACTOR, sep="_")
#'y <- na.omit(y)
#'n <- length(y$ID)
#'output <- CVrandomWithinFactor(y$ID, factorID=y$IDUnique, seed=123, 
#'k=5, exclusive=TRUE)
#'table(output)
#'head(output)
#'
CVrandomWithinFactor <- function(ID, factorID, seed, k, exclusive) {
  set.seed(seed)
  which.pop <- unique(factorID)
  y.u <- unique(ID)
  val.samp3 <- NULL
  for (j in 1:length(which.pop)) {
    y2 <- matrix(y.u[factorID == which.pop[j]], ncol = 1)
    set.seed(seed + j)
    modu <- nrow(y2) %% k
    if (!modu == 0) 
      val.samp <- sample(c(rep(1:k, each = (nrow(y2) - 
                                              modu)/k), sample(1:k, modu)), nrow(y2), replace = FALSE)
    if (modu == 0) 
      val.samp <- sample(rep(1:k, each = (nrow(y2))/k), 
                         nrow(y2), replace = FALSE)
    val.samp2 <- data.frame(y2, val.samp)
    val.samp3 <- as.data.frame(rbind(val.samp3, val.samp2))
  }
  scheme <- val.samp3[order(as.character(val.samp3[, 1])), "val.samp"]
  names(scheme) <- ID
  return(scheme)
}

#' Function to make a CV scheme based on a commited test set.
#'  
#' Function to make a CV scheme based on a commited test set.
#' 
#'@param trainingSet character vector of the observations in the training set.
#'@param validationSet character vector of the observations in the specified test set.
#'@param k integer value for the number of folds used in the k-cross-validation.
#'@return named vector of numeric scores showing the assignment of the observations 
#'to their respective set used in the k-fold cross-validation.
#'@export
#'@author Ruud Derijcker
#'@references Based on synbreed's crossVal function
#'
CVcommit <- function(trainingSet, validationSet, k) {
  replication <- length(names(trainingSet))
  k <- length(trainingSet[[1]])
  val.samp2 <- as.data.frame(y$ID)
  val.samp2$val.samp <- rep(NA, n)
  k <- length(names(trainingSet[[i]]))
  for (ii in 1:k) {
    val.samp2[val.samp2[, 1] %in% trainingSet[[i]][[ii]], 2] <- ii
  }
  val.samp3 <- val.samp2
  return(val.samp3)
}

#' Function to generate a table of combinations when applying permutation
#' 
#' Function to generate a table of combinations when applying permutation
#' 
#'@param n numeric for the set size where to apply the permutations
#'@param r numeric for the number of wished combinations for the permutations
#'@param v numeric value used to show the number of combinations in the output. 
#'The argument is passed to R's base permute function. We use the default of 1:n, 
#'which is the full set is outputted.
#'@return data frame with C(n,r) x r combinations.
#'@details we provide a table of the permutation combinations given by the number 
#'of r-combinations of an n-set, C(n,r).
#'@export
#'@examples
#'head(permute(10,2))
#'
permute <- function (n, r, v = 1:n) {
  if (r == 1) {
    matrix(v, n, 1)
  } else if (n == 1) {
    matrix(v, 1, r)
  } else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], permute(n -1, r - 1, v[-i])))
    X
  }
}

#' Function to make a CV scheme based on random sampling of entry IDs 
#' in an incomplete field trial setup.
#'  
#' Function to make a CV scheme based on random sampling of entry IDs 
#' in an incomplete field trial setup.
#' 
#'@param ID character vector of the observation IDs. 
#'The names are a combination of the entry name and the location.
#'@param factorID character vector of the entry IDs used in the randomization.
#'@param LOCATION character describing the column name representing the location 
#'information.
#'@param k integer value for the number of folds used in the k-cross-validation.
#'@param exclusive logical whether sampling should be done with replacement. The 
#'argument is passed to the replace argument of the samp.int function as the negation, 
#'i.e. exclusive is TRUE means replace=FALSE, such that the probability of 
#'choosing the next item is proportional to the weights amongst the remaining 
#'items. 
#'@param seed numeric value for the seed value used for the randomization by the 
#'set.seed function. In this way randomization can be reproduced by the user.
#'Default is NULL, which uses 123 as value for the seed.
#'@return named vector of numeric scores showing the assignment of the observations 
#'to their respective set used in the k-fold cross-validation.
#'@export
#'@details we developed our own functionality to simulate randomization for an 
#'incomplete field trial, which is based on permutation using the attached 
#'permute function. in the incomplete field trial setup we define an equal set 
#'of entries at every location, where we ensure no overlap in entry IDs 
#'during this process.
#'@author Ruud Derijcker
#'@examples
#'data(exampleCV)
#'y <- exampleCV[,which(colnames(exampleCV) %in% c("GERMPLASM", "LOCATION"))]
#'colnames(y) <- c("IDUnique","FACTOR")
#'y$ID <- paste(y$IDUnique, y$FACTOR, sep="_")
#'y <- na.omit(y)
#'n <- length(y$ID)
#'output <- CVincompleteTrial(ID=y$ID, factorID=y$IDUnique, LOCATION=y$FACTOR, 
#'                            seed=123, k=5, exclusive=TRUE)
#'table(output)
#'head(output)
#'
CVincompleteTrial <- function(ID, factorID, LOCATION, k, exclusive, seed) {
  set.seed(seed)
  locations <- as.character(unique(LOCATION))
  nloc <- length(locations)
  Levels <- levels(factor(factorID))
  nLevels <- length(Levels)
  n <- floor(nLevels / nloc)
  #modu <- length(ID) %% (k*nloc)
  #n <- ((length(ID) -  modu)/k)/nloc
  
  res <- array(dim=c(n, ncol=length(locations), k))
  res2 <- matrix(nrow=n, ncol=length(locations))
  df2 <- df3 <- NULL
  
  for(l in 1: length(locations)){
    plantNames <- as.character(factorID)[LOCATION==locations[l]]
    plantNames2 <- setdiff(plantNames, c(res2[, 1:l]))
    res2[,l] <- as.character(sample(plantNames2, n, replace=!exclusive))
  }
  tmp <- t(apply(permute(nloc,nloc), 1, function(x) diff(x)))
  select <- apply(tmp, 1, function(x) all(x %in% c(1,nloc-(2*nloc-1)))) #-4
  permScheme <- permute(nloc,nloc)[select,]
  for(i in 1: k)
    res[,,i] <- res2[,permScheme[i,]]
  
  for(i in 1: k) {
    for(l in 1: nloc){
      df2 <- data.frame(name=paste(res[,l,i],locations[l], sep="_"), fold=i)
      df3 <- rbind(df3,df2)
    }
  }
  names2 <- setdiff(ID, df3$name)
  if(length(names2)!=0) {
    remainders <- data.frame(name=setdiff(ID, df3$name), fold=k)
    results <- rbind(df3, remainders)
  } else {
    results <- df3
  }
  scheme <- results$fold
  names(scheme) <- results$name
  scheme <- scheme[match(ID, names(scheme))]
  return(scheme)
}

#' Cross-validation function
#' 
#' Cross-validation function used in combination with BGLR
#'@param x a data frame with at least the following information:\describe{
#' \item{\code{GERMPLASM}:}{Name of then entries.}
#' \item{\code{LOCATION}:}{Name of the geographic locations of the multi-field 
#' trial. In this function we assume the factor used in the setting up the 
#' cross-validation schemes is the geographic location. However this can be any 
#' field design factor which is adequate as analysis factor.}
#' }
#'@param id character specifying the column name of the entries IDs in x. 
#'Default is GERMPLASM.
#'@param factor character specifying the column name of the factor to use in the 
#'cross-validation in x. Default is LOCATION, refering to the graphical locations 
#'in considering a multi-location field trial.
#'@param k integer defining the number of folds for k-fold cross validation, 
#'thus k should be in [2,nrow(y)], where y is the vector of phenotypic values. 
#'The default is 5.
#'@param replication numeric defining the number of replications of the 
#'cross-validation. Default is 3. 
#'@param seed numeric value for the seed value used for the randomization by the 
#'set.seed function. In this way randomization can be reproduced by the user. 
#'Default is NULL, which uses 123 as value for the seed.
#'@param exclusive logical whether sampling should be done with replacement. The 
#'argument is passed to the replace argument of the samp.int function as the negation, 
#'i.e. exclusive is TRUE means replace=FALSE, such that the probability of 
#'choosing the next item is proportional to the weights amongst the remaining 
#'items.
#'@param sampling character specifying which sampling strategy to use in the 
#'cross-validation. The different sampling strategies are described below:
#'\describe{
#'      \item{\code{randomByID}:}{Random sampling by name of the observations}
#'      \item{\code{randomAccrossFactor}:}{Random sampling by name of the entries 
#'      taking into account randomization across a defined factor}
#'      \item{\code{randomByFactor}:}{Random sampling by name of the entries 
#'      using the factor to define the sets}
#'      \item{\code{randomWithinFactor}:}{Random sampling by name of the entries 
#'      taking into account randomization within a defined factor}
#'      \item{\code{popStructureAccrossFactor}:}{Accounts for across population 
#'      structure information, e.g. test and training sets contain a set of 
#'      complete families}
#'      \item{\code{popStructureWithinFactor}:}{Accounts for within population 
#'      structure information, e.g. each family is splitted into k subsets}
#'      \item{\code{commit}:}{Sampling done using defined test and training sets}
#'      \item{\code{incompleteTrial}:}{Random sampling by taking into account 
#'      an incomplete field trial setup}
#'}
#'If sampling is "commit" the sets of names have to specified in the trainingSet 
#'and validationSet arguments.
#'@param trainingSet character vector of the observations in the training set.
#'@param validationSet character vector of the observations in the specified test 
#'set.
#'@param populationStructure vector of length nrow(y) assigning individuals to 
#'a population structure, where y refers to the vector of phenotypes. This 
#'argument is only required for the options sampling="popStructureAccrossFactor" 
#'or sampling="popStructureWithinFactor".
#'@param verbose logical whether to output information about the progress of the 
#'cross-validation. Default is FALSE.
#'@details in cross validation (CV) the data set is splitted into a training set, 
#'and a validation or test set. For sampling into the sets, k-fold cross validation 
#'is applied, where the data set is splitted into k subsets and k-1 comprising 
#'the training set and 1 is the test set, repeated for each subset. The function 
#'is based on the crossVal function from the synbreed package. We made the 
#'function more flexible by taking out the cross-validation schemes functionality, 
#'to allow easy plug-in of more user-defined CV schemes. Further, the function 
#'was adjusted to work with the BGLR framework.
#'@return data frame with the result of the sampling of the entries into k-folds 
#'using a number of user-defined replications. The table includes following columns:
#'\itemize{
#'  \item{ID}{The names of the observations.}
#'  \item{Rep[x]}{[x] columns of numeric scores according to the assignment of the 
#'  observations into 1...k folds, where [x] is set by the replication argument}
#'}
#'@references \describe{
#'      \item{\code{1}:}{Albrecht T, Wimmer V, Auinger HJ, Erbe M, Knaak C, 
#'      Ouzunova M, Simianer H, Schoen CC (2011) Genome-based prediction of 
#'      testcross values in maize. Theor Appl Genet 123:339-350.}
#'      \item{\code{2}:}{Gustavo de los Campos and Paulino Perez Rodriguez 
#'      (2014). BGLR: Bayesian Generalized Linear Regression. R package version 
#'      1.0.3. http://CRAN.R-project.org/package=BGLR}
#'}
#'@export
#'@examples
#' data(exampleCV)
#' scheme1 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomByFactor",verbose=TRUE)
#' scheme2 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION",
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="incompleteTrial",verbose=TRUE)
#' scheme3 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomAccrossFactor",verbose=TRUE)
#' scheme4 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION",
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomWithinFactor",verbose=TRUE)  
#' scheme5 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomByID",verbose=TRUE)
#' head(scheme1)
#' head(scheme2)
#' head(scheme3)
#' head(scheme4)
#' head(scheme5)
#' 
crossValidate <- function(x, id="GERMPLASM", factor="LOCATION",
                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
                          sampling=c("randomByID","randomAccrossFactor","randomByFactor",
                                     "randomWithinFactor", "popStructureAccrossFactor",
                                     "popStructureWithinFactor", "commit", "incompleteTrial"),
                          trainingSet=NULL, validationSet=NULL, 
                          populationStructure=NULL, verbose=FALSE) {
  
  # prepare to take into account a factor for sampling
  if(!is.null(factor)) {
    y <- x[,which(colnames(x) %in% c(id, factor))]
    colnames(y) <- c("IDUnique","FACTOR")
    y$ID <- paste(y$IDUnique, y$FACTOR, sep="_")
  } else {
    y <- x[,which(colnames(x) %in% c(id))]
    colnames(y) <- c("ID")
  }
  y <- na.omit(y)
  n <- length(y$ID)
  
  # some checks
  if (k < 2) 
    stop("folds should be equal or greater than 2")
  if (k > n) 
    stop("folds should be equal or less than the number of observations")
  
  # sampling schemes
  sampling <- match.arg(sampling)
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    seed <- 123
    set.seed(seed)
  }
  
  if(sampling == "randomByID") {
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomByID(ID=y$ID, k=k, exclusive=exclusive, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if (sampling == "randomAccrossFactor" & !is.null(factor)) {
    if (any(is.na(y$FACTOR))) 
      stop("no missing values allowed in factor")
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomAccrossFactor(ID=y$ID, factorID=y$IDUnique, exclusive=exclusive, k=k, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if (sampling == "randomByFactor") {
    if (any(is.na(y$FACTOR))) 
      stop("no missing values allowed in factor")
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomAccrossFactor(ID=y$ID, factorID=y$FACTOR, exclusive=exclusive, k=k, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if (sampling == "randomWithinFactor") {
    if (any(is.na(y$FACTOR))) 
      stop("no missing values allowed in factor")
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomWithinFactor(ID=y$ID, factorID=y$IDUnique, exclusive=exclusive, k=k, seed=x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if (sampling == "incompleteTrial") { 
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVincompleteTrial(ID=y$ID, factorID=y$IDUnique, LOCATION=y$FACTOR, 
                        exclusive=exclusive, k=k, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if(sampling == "popStructureAccrossFactor") {
    if (length(populationStructure) != length(x$ID)) 
      stop("population structure must have equal length as obsersvations in data")
    if (any(is.na(y$FACTOR))) 
      stop("no missing values allowed in factor")
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomAccrossFactor(ID=y$ID, factorID=y$populationStructure, 
                            exclusive=exclusive, k=k, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if(sampling == "popStructureWithinFactor") {   
    if (any(is.na(y$FACTOR))) 
      stop("no missing values allowed in factor")
    scheme <- data.frame(ID=y$ID, apply(data.frame(1:replication), 1, function(x) 
      CVrandomWithinFactor(ID=y$ID, factorID=y$populationStructure, 
                           exclusive=exclusive, k=k, seed=seed+x)))
    colnames(scheme) <- c("ID", paste("Rep", 1:replication, sep=""))
    if(verbose)
      data.frame(Fold=paste("fold",1:k, sep=""),apply(scheme[,-1], 2, table))
    
  } else if (sampling == "commit") {
    if(is.null(trainingset))
      stop("trainingSet have to be specified")
    stop("under construction")
    scheme <- CVcommit (trainingSet, validationSet, k)
  }
  return(scheme)
}

#' Plot a CV scheme based on crossValidate function
#' 
#' Function to graphically show a CV scheme made by the crossValidate function
#'@param data a data frame with at least the following information:\describe{
#' \item{\code{GERMPLASM}:}{Name of then entries.}
#' \item{\code{LOCATION}:}{Name of the geographic locations of the multi-field 
#' trial.}
#' }
#'@param scheme data frame output from the crossValidate function which was based 
#'on a user decided sampling strategy to use in the cross-validation. 
#'@param title character defining the title of the plot
#'@return plot based on the ggplot package which shows the splitting of the training 
#'and validation sets for the first replication. Further, the plot is split by the 
#'levels of the chosen factor, which is generally the location. 
#'@export
#'@author Ruud Derijcker
#'@examples
#'data(exampleCV)
#' scheme1 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomByFactor",verbose=TRUE)
#' scheme2 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION",
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="incompleteTrial",verbose=TRUE)
#' scheme3 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomAccrossFactor",verbose=TRUE)
#' scheme4 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION",
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomWithinFactor",verbose=TRUE)  
#' scheme5 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomByID",verbose=TRUE)
#'
#'plotCV(exampleCV, scheme1, "CV-scheme based on randomByFactor")
#'plotCV(exampleCV, scheme2, "CV-scheme based on incompleteTrial")
#'plotCV(exampleCV, scheme3, "CV-scheme based on randomAccrossFactor")
#'plotCV(exampleCV, scheme4, "CV-scheme based on randomWithinFactor")
#'plotCV(exampleCV, scheme5, "CV-scheme based on randomByID")
#'
plotCV <- function(data, scheme, title){
  data <- data.frame(data, scheme)
  data <- data[order(data$ROW, data$RANGE),]
  data$"FoldScheme" <- as.factor(ifelse(data$Rep1==1,"training", "validation"))
  pg1 <- ggplot(data, aes(ROW, RANGE, fill = FoldScheme)) +
    facet_grid(~ LOCATION) + geom_tile() + labs(title=title)
  plot(pg1)
}

#' Add observations to a data frame so that the design is rectangular (per
#' field)
#' 
#' Add observations to a data frame so that the design is rectangular (per
#' field).
#' 
#' For each level of \code{fields}, a rectangular design is created by adding
#' empty observations for the missing \code{rows} and \code{ranges}.
#' 
#' @param object A data frame with at least the columns \code{rows}, \code{ranges}
#' and \code{fields}, where fields refer to locations.
#' @param fields Character string specifying the column name with the fields names 
#' corresponding to the field locations in \code{object}.
#' @param rows (Optional) Character string specifying the column name with the
#' ROW-coordinates in \code{object}.  Default = \code{"ROW"}.
#' @param ranges (Optional) Character string specifying the column name with
#' the RANGE-coordinates in \code{object}.  Default = \code{"RANGE"}.
#' @return An extended version of the data frame \code{object}.
#' @author Ruud Derijcker
#' @export
#' @examples
#' data(exampleCV)
#' scheme1 <- crossValidate(x=exampleCV, id="GERMPLASM", factor="LOCATION", 
#'                          k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                          sampling="randomByFactor",verbose=TRUE)
#' P2 <- data.frame(exampleCV, scheme1)
#' P2 <- make_design_rectangular(P2, fields="LOCATION",rows="ROW",ranges="RANGE")
#' head(P2)
#' 
make_design_rectangular <- function(object, fields, rows="ROW", ranges="RANGE") {
  if(any(!c(fields, rows, ranges) %in% colnames(object)))
    stop(paste("The column(s)", paste(setdiff(c(fields, rows, ranges), 
                                              colnames(object)), collapse = ", "), 
               "is/are not found in 'x'."))
  
  object$tmp <- factor(object[, fields])
  object$tmp1 <- object[, rows]
  object$tmp2 <- object[, ranges]
  
  out <- NULL
  for(i in levels(object$tmp)){
    sub <- subset(object, tmp == i)
    min.range <- min(sub$tmp2)
    max.range <- max(sub$tmp2)
    min.row <- min(sub$tmp1)
    max.row <- max(sub$tmp1)
    
    full.design <- expand.grid(tmp2 = min.range:max.range, tmp1 = min.row:max.row)
    missing <- subset(merge(sub, full.design, by = c("tmp2", "tmp1"), all.y = TRUE),
                      is.na(tmp))
    if(nrow(missing) > 0) {
      sub$tmp <- as.character(sub$tmp)
      missing$tmp <- as.character(missing$tmp)
      missing$tmp <- i
      out <- rbind(out, sub, missing)
    }
    else {
      sub$tmp <- as.character(sub$tmp)
      out <- rbind(out, sub)
    }
  }
  out[, rows] <- out$tmp1
  out[, ranges] <- out$tmp2
  out[, fields] <- factor(out$tmp)
  return(out[, 1:(ncol(object) - 3)])
}

#' Plot a CV scheme for the data used in the thesis
#' 
#' Plot a CV scheme for the data used in the thesis
#'@param data a data frame with at least the following information:\describe{
#' \item{\code{GERMPLASM}:}{Name of then entries.}
#' \item{\code{LOCATION}:}{Name of the geographic locations of the multi-field 
#' trial.}
#' \item{\code{ROW}:}{Name of then rows}
#' \item{\code{RANGE}:}{Name of then ranges}
#' }
#'@param scheme data frame output from the crossValidate function which was based 
#'on a user decided sampling strategy to use in the cross-validation. 
#'@param title character defining the title of the plot
#'@return plot based on the ggplot package which shows the splitting of the training 
#'and validation sets for the first replication. Further, the plot is split by the 
#'levels of the chosen factor, which is generally the location. 
#'@export
#'@author Ruud Derijcker
#'@examples
#'data(P)
#'scheme1 <- crossValidate(x=P, id="GERMPLASM", factor="LOCATION", 
#'                         k=5, replication=3, seed=NULL, exclusive=TRUE, 
#'                         sampling="incompleteTrial",verbose=TRUE)
#'plotCVThesisData(P, scheme1, "CV scheme 'incompleteTrial' on thesis data set")
#'
plotCVThesisData <- function(data, scheme, title="CV scheme") {
  P2 <- data.frame(data, scheme)
  P2 <- make_design_rectangular(P2, fields="LOCATION",rows="ROW",ranges="RANGE")
  P2 <- P2[order(P2$LOCATION, P2$GERMPLASM),]
  P2$CVrow <- rep(1:10, each=33)
  P2$CVrange <- rep(1:33, each=1)
  P2$cat <- ifelse(P2$Rep1==1,1,0)
  P2 <- P2[order(P2$GERMPLASM, P2$CVrow, P2$CVrange),]
  pg1 <- ggplot(P2, aes(CVrow, CVrange, fill = factor(cat))) +
    facet_grid(~ LOCATION) + geom_tile() + labs(title=title)
  plot(pg1)
}
