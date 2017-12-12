#' Create a field-layout plot from a data frame.
#' 
#' Create a field-layout plot from a data frame.
#' 
#' @param x A data frame containing at least the columns \code{column},
#' \code{rows} and \code{ranges}.
#' @param column Character string specifying the name of the column in \code{x}
#' to be used for coloring the different positions in the field.  This should
#' be either a factor-column or a numeric column.
#' @param textColumn (Optional) Character string specifying the name of the column in 
#' \code{x} to be used to print text within the plots
#' @param asFactor Logical indicating whether the specified column should be 
#' considered as a factor, even if it is numeric.
#' @param borderColumns (Optional) A chracter vector specifying the factors 
#' (maximum of three factors) to be used to draw border columns in the field map. 
#' Default = \code{c("EXPT", "BLOCK")}.
#' The factors are assumed to be given in nested order, such as \code{c("EXPT", "TRIAL", "BLOCK")}, 
#' \code{c("EXPT", "BLOCK", "SUBBLOCK")} etc. 
#' -- if only one factor is given, a bold border will be drawn
#' -- if two factors are given, a bold border will be drawn for the main factor with regular border 
#' for the sub-factor
#' -- if three factors are given, a bold border will be drawn for the main factor, regular border 
#' will be drawn for the sub-factor and dotted-line border will be drawn for sub-sub-factor
#' @param colors Charachter vector speciyfing the colors to be used for the
#' different levels in \code{column}.  Default = \code{NULL} in which case
#' pre-specified colors are used.
#' @param colorscale (Optional) A numeric vector of length 2 with the range to
#' use the \code{colors} over and is applicable only when 'column' is numeric variable. 
#' Default is \code{NULL} in which case the observed range is computed from the \code{column}.
#' @param rows Character string specifying the column name with the
#' ROW-coordinates in \code{x}.
#' @param ranges Character string specifying the column name with the
#' RANGE-coordinates in \code{x}.
#' @param \dots Extra graphical parameters for the
#' \code{\link[lattice]{levelplot}}-function.
#' @return Plot is submitted to the graphical device.
#' @author Ruud Derijcker, Katrien Baert, Bert Lataire
#' @export
#' @examples
#' data(exampleCV)
#' exampleCV$CHECK <- sample(c(rep(0, 80), 
#' sample(c(1:4), nrow(exampleCV), 20)), 100)
#' exampleCV$EXPT <- exampleCV$BLOCK <- NA
#' plotFields(x=exampleCV, column="CHECK", main="Fieldplot for dummy example",
#'            textColumn="GERMPLASM",
#'            colors=c("grey", "red","blue","green","yellow"))
#'
plotFields <- function(x, column, textColumn = NULL, asFactor=FALSE,
                       borderColumns=c("EXPT", "BLOCK"), colorscale=NULL,
                       colors = NULL, rows = "ROW", ranges = "RANGE", ...) {
  require(lattice)
  if(length(borderColumns) > 3)
    stop("The 'borderColumn' argument can only take up to three factors. \n", 
         "Please check the argument.")
  if(!all(borderColumns %in% colnames(x))) {
    missingColumns <- setdiff(borderColumns, colnames(x))
    warning("The following borders column(s) is/are not found in 'x': ", 
            paste(missingColumns, collapse = ", "))
    borderColumns <- borderColumns[-which(borderColumns %in% missingColumns)]
  }
  if(!rows %in% colnames(x))
    stop("The row-indicator 'rows' is not found as column of 'x'.")
  if(!ranges %in% colnames(x))
    stop("The range-indicator 'ranges' is not found as column of 'x'.")
  x$ROW <- x[, rows]
  if(!is.numeric(x$ROW))
    x$ROW <- as.numeric(as.character(x$ROW))
  x$RANGE <- x[, ranges]
  if(!is.numeric(x$RANGE))
    x$RANGE <- as.numeric(as.character(x$RANGE))
  
  x$tmpColumn <- "FIELD1"
  x <- make_design_rectangular(x, fields = "tmpColumn")
  data <- x
  data$Y <- data[, column]
  if(asFactor)
    data$Y <- factor(data$Y)
  
  if(is.numeric(data$Y)){
    if(is.null(colors))
      colors <- rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", 
                      "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", 
                      "#74ADD1", "#4575B4", "#313695"))
    if(!is.null(colorscale)){
      breakpoints <- seq(colorscale[1], colorscale[2], length = 11)
    } else {
      if(all(c(-1, 1) %in% sign(data$Y))) {
        M <- max(abs(data$Y), na.rm = TRUE)
        breakpoints <- seq(-M, M, length = 11)
      } else {
        breakpoints <- seq(min(data$Y, na.rm = TRUE), max(data$Y, na.rm = TRUE), length = 11)
      }
    }  	
    colorkey <- TRUE
    border <- "white"
  } else {
    levels <- levels(factor(data$Y))
    if(is.null(colors))
      colors <- colorRampPalette(c("cyan","green", "blue", "red", "yellow"))
    breakpoints <- 0:length(levels)
    colorkey <- list(labels = list(labels = levels, at = 0:length(levels) + 0.5))
    border <- "white"
  }
  
  if(!is.null(textColumn)){
    data$textColumn <- data[, textColumn]
    data$textColumn <- ifelse(is.na(data$textColumn), "", as.character(data$textColumn))
  }
  
  nBorderColumns <- length(borderColumns)
  if(nBorderColumns >= 1) {
    data$border1 <- factor(data[, borderColumns[1]])
  }
  if(nBorderColumns >= 2) {
    data$border2 <- factor(data[, borderColumns[2]])
  }
  if(nBorderColumns == 3) {
    data$border3 <- factor(data[, borderColumns[3]])
  }
  
  # Actual plotting
  require(lattice)
  levelplot(Y ~ ROW*RANGE, 
            data = data, 
            contour = FALSE, 
            colorkey = colorkey,
            at = breakpoints,
            col.regions = colors, 
            border = border, 
            panel = function(x, y, z, ...){
              panel.levelplot(x, y, z,...)
              if(!is.null(textColumn)){
                panel.text(x = x, y= y, data$textColumn)
              }
              if("border1" %in% colnames(data)){
                for(iBorder1 in levels(data$border1)){
                  sub1 <- subset(data, border1 == iBorder1)
                  panel.rect(min(sub1$ROW - 0.5), min(sub1$RANGE - 0.5), 
                             max(sub1$ROW + 0.5), max(sub1$RANGE + 0.5), 
                             lwd = 2, lty = 1, border = "black")
                }          
              }
              if(all(c("border1", "border2") %in% colnames(data))){
                for(iBorder1 in levels(data$border1)){
                  for(jBorder2 in levels(data$border2)){
                    sub2 <- subset(data, border1 == iBorder1 & border2 == jBorder2)
                    if(nrow(sub2) > 0){
                      panel.rect(min(sub2$ROW - 0.5), min(sub2$RANGE - 0.5), 
                                 max(sub2$ROW + 0.5), max(sub2$RANGE + 0.5), 
                                 lwd = 1, lty = 1, border = "black")
                    }              
                  } 
                }				           
              }
              if(all(c("border1", "border2", "border3") %in% colnames(data))){
                for(iBorder1 in levels(data$border1)){
                  for(jBorder2 in levels(data$border2)){
                    for(kBorder3 in levels(data$border3)){
                      sub3 <- subset(data, border1 == iBorder1 & border2 == jBorder2 
                                     & border3 == kBorder3)
                      if(nrow(sub3) > 0){
                        panel.rect(min(sub3$ROW - 0.5), min(sub3$RANGE - 0.5), 
                                   max(sub3$ROW + 0.5), max(sub3$RANGE + 0.5), 
                                   lwd = 1, lty = 2, border = "black")
                      }				      
                    }              
                  } 
                }				           
              }
              
            }, ...)		
}

#' Boxplots with checks alongside the boxes
#' 
#' Create boxplots for a sequence of traits for each level of a grouping
#' factor.
#' 
#' @param x A data frame containing at least the columns \code{trait} and
#' \code{factor} and also \code{checks} and \code{checks_factor} if both
#' specified.
#' @param traits Charachter vector with the column names in \code{data} of the
#' traits to be plotted.
#' @param factor Character string with the column name in \code{data} of the
#' factor by which to plot the trait.
#' @param checks (Optional) Charachter vector with the names of the checks in
#' the \code{checks_factor} column.
#' @param checks_factor (Optional) Charachter string with the column name in
#' \code{data} of the factor that contains the observations to be plotted as
#' \code{checks}.
#' @param plottype (Optional) Character string speciyfing which type of plot to
#' produce, should be either \code{"Boxplot"} or \code{"Stripplot"}. Default =
#' \code{"Boxplot"}.
#' @param boxcolor (Optional) Character vector specifying the color(s) for the
#' boxplots. Default = \code{"white"}.
#' @param stripcolor (Optional) Character vector specifying the color(s) for
#' the points in the stripplot. Default = \code{"red"}.
#' @param checks.color (Optional) Character string specifying the color(s) for
#' the checks. Default = \code{"blue"}.
#' @param legend (Optional) Logical indicating whether a legend for the
#' \code{checks} should be added to the plot. Default = \code{TRUE}.
#' @param trim.labels (Optional) Numeric value specifying how many characters
#' of the levels of \code{factor} should be plotted on the X-axis.  If
#' \code{NULL}, the levels are plotted as they are in \code{x}. Default =
#' \code{NULL}.
#' @param \dots Extra graphical parameters for the plot-function.
#' @return Plot is submitted to the graphical device.
#' @author Ruud Derijcker, Katrien Baert
#' @export
#' @examples
#' 

#' data(exampleCV)
#' exampleCV$CHECK <- sample(c(rep(0, 300), 
#' sample(c(1:4), nrow(exampleCV), 20)), 100)
#' exampleCV$EXPT <- exampleCV$BLOCK <- NA
#' colnames(exampleCV)[colnames(exampleCV)== "trait"] <- "dummy"
#' exampleCV$dummy <- abs(exampleCV$dummy)*100
#' trait <- "dummy"
#' boxplots_with_checks(x=exampleCV, traits="dummy", factor="LOCATION", 
#'       checks=exampleCV$GERMPLASM[exampleCV$CHECK!=0], 
#'       checks_factor = "GERMPLASM", checks.color = "grey", 
#'      boxcolor="blue", main="Yield by location\n(all checks)")
#'
boxplots_with_checks <- function(x, traits, factor, checks=NULL, checks_factor=NULL,
                                 plottype="Boxplot", boxcolor="white",
                                 stripcolor="red", checks.color="blue",
                                 legend=TRUE, trim.labels=NULL, ...) {
  # Checks
  if(!all(traits %in% colnames(x)))
    stop("Not all trait are found in 'x'.")
  
  if(!factor %in% colnames(x))
    stop("The specified 'factor' is not found in the 'data' frame.")
  
  x[, factor] <- factor(x[, factor])
  levels <- levels(x[, factor])
  n.levels <- length(levels)
  
  if(!is.null(checks_factor) & !is.null(checks)){
    if(!checks_factor %in% colnames(x))
      stop("The specified 'checks_factor' is not found in the 'data' frame.")
    
    if(!all(checks %in% x[, checks_factor]))
      stop("Not all 'checks' are levels of the 'check_factor'.")
    
    x[, checks_factor] <- factor(x[, checks_factor])
    sub <- x[!x[, checks_factor] %in% checks, ]
    sub_checks <- x[x[, checks_factor] %in% checks, ]
  } else {
    checks <- NULL
    sub <- x
  }
  
  # Plot
  for(i in 1:length(traits)){
    if(!trait %in% colnames(x))
      stop(paste("The trait", trait, "is not found in 'x'."))
    
    if(!factor %in% colnames(x))
      stop("The specified 'factor' is not found in the 'data' frame.")
    
    x[, factor] <- factor(x[, factor])
    levels <- levels(x[, factor])
    n.levels <- length(levels)
    
    if(!is.null(checks_factor) & !is.null(checks)){
      if(!checks_factor %in% colnames(x))
        stop("The specified 'checks_factor' is not found in the 'data' frame.")
      
      if(!all(checks %in% x[, checks_factor]))
        stop("Not all 'checks' are levels of the 'check_factor'.")
      
      x[, checks_factor] <- factor(x[, checks_factor])
      sub <- x[!x[, checks_factor] %in% checks, ]
      sub_checks <- x[x[, checks_factor] %in% checks, ]
    } else {
      checks <- NULL
      sub <- x
    }
    
    dots <- as.list(substitute(list(...)))[-1]
    
    # Plot
    formula <- as.formula(paste(trait, "~", factor))
    if(plottype == "Boxplot"){
      b <- boxplot(formula, data = sub, col = boxcolor, xaxt = "n", ...)
    } else if(plottype == "Stripplot"){
      b <- boxplot(formula, data = sub, plot = FALSE, ...)
      stripchart(formula, data = sub, col = stripcolor, 
                 method = "jitter", vertical = T, ...)
    }
    
    # Add (trimmed) labels to the x-axis
    if(is.null(trim.labels))
      axis(side = 1, at = 1:n.levels, labels = levels, ...)
    else
      axis(side = 1, at = 1:n.levels, labels = substr(levels, 1, trim.labels), ...)
    
    # Add default title and y-label
    if(is.null(eval(dots$main)))
      title(main = trait)
    if(is.null(eval(dots$ylab)))
      title(ylab = trait)
    
    # Add checks
    if(!is.null(checks)) {
      points(as.numeric(sub_checks[, factor]) + 0.2, sub_checks[, trait], 
             col = checks.color, pch = 19)
      
      if(legend){
        par(xpd = TRUE)
        usr <- par("usr")
        legend(usr[1], usr[4] + diff(usr[3:4])*0.07, 
               "Checks", pch = 19, col = checks.color, bty = "n")
        par(xpd = FALSE)
      }
    }
  }
}

#' Likelihood ratio test for two asreml-models.
#' 
#' Perform a likelihood ratio test for two asreml-models that have an equal
#' fixed part and a nested random part.
#' 
#' The likelihood ratio test is performed by comparing two times the difference
#' in log-likelihood with the quantiles if the \eqn{\chi^2}-distribution with
#' as degrees of freedom the difference in number of variance components
#' estimated in both models.
#' 
#' If the degrees of freedom equal 1, the obtained p-value is divided by 2.
#' 
#' @param model0 An asreml-object, see \code{\link[asreml]{asreml}}.
#' @param model1 An asreml-object with the same fixed part as \code{model0} and
#' a random part that is an extension of the random part of \code{model0}.
#' @return A list with class "\code{htest}" containing the following
#' components: \describe{ \item{statistic}{The value of the test statistic.}
#' \item{parameter}{The degrees of freedom for the test statistic.}
#' \item{p.value}{The p-value for the test.} \item{conf.int}{\code{NULL}}
#' \item{estimate}{\code{NULL}} \item{null.value}{\code{NULL}}
#' \item{alternative}{\code{Extended model outperforms the 'null model'.}}
#' \item{method}{\code{Likelihood ratio test (Chi-square)}}
#' \item{data.name}{Names of the two models that were tested.} }
#' @author Ruud Derijcker, Katrien Baert
#' @export 
#' @examples
#' 
#' library(asreml)
#' data(oats)
#' model0 <- asreml(fixed = yield ~ Variety*Nitrogen, random = ~ Blocks:Wplots, 
#'           data = oats)
#' model1 <- asreml(fixed = yield ~ Variety*Nitrogen, random = ~ Blocks/Wplots, 
#'           data = oats)
#' LR.test.asreml(model0, model1)
#' 
LR.test.asreml <- function(model0, model1) {
  # Check the class of the input arguments
  if(!(class(model0) == "asreml" && class(model1) == "asreml"))
    stop("The input arguments should be asreml-objects.")
  # Check equality of fixed effects
  if(model0$fixed.formula != model1$fixed.formula)
    stop("Both models should have the same fixed part.")
  
  v0 <- summary(model0)$varcomp
  v0 <- v0[v0$constraint != "Fixed", ]
  v1 <- summary(model1)$varcomp
  v1 <- v1[v1$constraint != "Fixed", ]
  
  # Check that the models have either a different random component OR error structure 
  if((nrow(v0) == nrow(v1)) && 
       (length(model0$R.param) == length(model1$R.param)))
    stop("Both models should have a different (nested!) random part or error 
         structure.")
  
  if(nrow(v0) > nrow(v1)){
    temp <- model0
    model0 <- model1
    model1 <- temp
    
    v0 <- summary(model0)$varcomp
    v0 <- v0[v0$constraint != "Fixed", ]
    v1 <- summary(model1)$varcomp
    v1 <- v1[v1$constraint != "Fixed", ]
  }
  
  # Check that random parts are nested
  if(any(!setdiff(rownames(v0), "R!variance") %in% rownames(v1))) 
    stop("model0 should have nested random part of model1.")
  
  df0 <- nrow(v0)
  df1 <- nrow(v1)	
  ll0 <- model0$loglik
  ll1 <- model1$loglik
  
  tstat <- 2*(ll1 - ll0)
  names(tstat) <- "X"
  df <- df1 - df0
  names(df) <- "degrees of freedom"
  p.val <- 1 - pchisq(tstat, df = df)
  if(df == 1)
    p.val <- p.val/2
  
  out <- list(statistic = tstat, parameter = df, p.value = p.val, 
              conf.int = NULL, estimate = NULL, null.value = NULL, 
              alternative = "Extended model outperforms the 'null model'.", 
              method = "Likelihood ratio test (Chi-square)", 
              data.name = paste(deparse(substitute(model0)), "and", 
                                deparse(substitute(model1))))
  class(out) <- "htest"
  return(out)
}

#' Plot variogram from residuals
#' 
#' @param x A data frame that contains at least the columns \code{residuals},
#' \code{ROW} and \code{RANGE}.
#' @param trim Logical indicating whether to trim the variogram to shorter lag 
#' values or not. Default = \code{FALSE}.
#' @param title Character string specifying the title for the plot. 
#' Default = \code{"Variogram"}.
#' @return Plot is submitted to the graphical device.
#' @author Ruud Derijcker, Krishna Bondalapati
#' @export 
plotVariogram <- function(x, trim=FALSE, title="Variogram"){ 
  
  if(!"package:asreml" %in% search()) {
    require(asreml)
  }
  vv <- asreml.variogram(x[, c("ROW", "RANGE", "Residuals")])
  if(!"x" %in% colnames(vv)) # Field is not nicely rectangular
    return()
  if(trim)
    vv <- vv[vv$np >= median(vv$np), ]
  
  p1 <- wireframe(gamma ~ x * y, 
                  data=vv, drape=F, colorkey=FALSE, screen=list(z=30, x=-60), 
                  aspect=c(1, 0.66), 
                  scales=list(distance=c(1.2, 1.2, 0.5), arrows=F, cex=0.65, col="black"), 
                  zlab="", xlab=list(label=paste("ROWf(lag)"), cex=0.65), 
                  ylab=list(label=paste("RANGEf(lag)"), cex=0.65), main=title,
                  par.settings=list(axis.line=list(col="transparent")))
  return(p1)
}

#' Plot a series of residual plots obtained from a linear mixed model
#' 
#' @param x A data frame that contains at least the columns \code{residuals},
#' \code{ROW}, \code{RANGE} and \code{by}.
#' @param residuals Character string specifying the column name of the
#' residuals in \code{x}.
#' @param by Character string specifying the column name for the splitting
#' factor in \code{x}, usually this is a location factor.
#' @param rows Character string specifying the column name with the \code{ROW}-
#' coordinates in x. Default = "ROW".
#' @param ranges Character string specifying the column name with the 
#' \code{RANGE}-coordinates in x. Default = "RANGE".
#' @param trimVariogram Logical indicating whether to trim the variogram to 
#' shorter lag values. Default = "FALSE".
#' @param trait Character string specifying the trait name to use in the title for
#' the plots. Default to \code{NULL}, in which case no secondary title with trait name
#' is printed. 
#' @return Produces various diagnostic residual plots. The colors are scaled to be 
#' comparable over different \code{by}-levels (if there are). Four plots are 
#' outputted per location and are shown in the graphical device. The plots are 
#' listed below and detailed about the construction and use is given in 
#' Appendix A in the referred thesis. \describe{ 
#'  \item{plot upper-left}{Variogram of the residuals for the rows and ranges. 
#'  By default the lags are limited by 30. }
#'  \item{plot upper-right}{Conditional residual plots, colored by their value and 
#'  sign.}
#'  \item{plot bottom-left}{Scatterplot of the the residuals against the range 
#'   number conditional on the row number}
#'  \item{plot bottom-right}{Scatterplot of the the residuals against the row 
#'   number conditional on the range number}
#' } 
#' @author Ruud Derijcker, Krishna Bondalapati
#' @references Derijcker, R. (2015). Investigating incorporation of 
#' genotype x environment interaction (G x E) for genomic selection in a 
#' practical setting. Unpublished M.Sc. thesis. University of Ghent: Belgium.
#' @export 
#' @seealso \code{\link{plotVariogram}}.
#' 
plotResiduals <- function(x, residuals, by, rows="ROW", ranges="RANGE", 
                          trimVariogram=FALSE, trait=NULL) {
  
  # CHECKS
  if(!residuals %in% colnames(x))
    stop("The column 'residuals' is not found in 'x'.")
  x$Residuals <- x[, residuals]
  if(!is.null(by) && !by %in% colnames(x))
    stop("The column 'by' is not found in 'x'.")
  if(any(!c(rows, ranges) %in% colnames(x)))
    stop(paste("The column(s)", paste(setdiff(c(rows, ranges), colnames(x)), 
                                      collapse=", "), "is/are not found in 'x'."))
  x$ROW <- x[, rows]
  x$RANGE <- x[, ranges]
  
  # PREPARE
  x$Group <- factor(x[, by])
  locs <- levels(x$Group)
  my.lim0 <- max(abs(x$Residuals), na.rm=TRUE)
  my.lim <- 1.1*my.lim0
  
  for(iLocation in locs) {
    sub <- subset(x, Group == iLocation)
    if(is.null(trait)) {
      title <- paste("Residuals at ", by, ": ", iLocation, sep="")  
    } else {
      title1 <- paste("Residuals at ", by, ": ", iLocation, sep="")
      title2 <- paste("TRAIT", trait, sep=": ")
      title <- paste(title1, "\n", title2)
    }
    print(plotVariogram(x=sub, trim=trimVariogram, title=title), 
          more=TRUE, position=c(xmin=0, ymin=0.5, xmax=0.5, ymax=1))
    print(plotFields(x=sub, column="Residuals",
                     main=title, colorscale=c(-my.lim0, my.lim0)),
          more=TRUE, position=c(xmin=0.5, ymin=0.5, xmax=1, ymax=1))
    print(xyplot(Residuals ~ RANGE, data=sub,  main=title, 
                 ylim=c(-my.lim0, my.lim0),
                 panel=function(x, y, ...){
                   panel.loess(x[!is.na(y)], y[!is.na(y)], col="red")
                   panel.abline(h=0, col="grey", lty=1)
                   panel.abline(lm(y ~ x), col="blue", lty=2)
                   panel.xyplot(x, y, col="black", ...)
                 }), more=TRUE, position=c(xmin=0, ymin=0, xmax=0.5, ymax=0.5))
    print(xyplot(Residuals ~ ROW, data=sub, main=title,
                 ylim=c(-my.lim0, my.lim0),
                 panel=function(x, y, ...){
                   panel.loess(x[!is.na(y)], y[!is.na(y)], col="red")
                   panel.abline(h=0, col="grey", lty=1)
                   panel.abline(lm(y ~ x), col="blue", lty=2)
                   panel.xyplot(x, y, col="black", ...)
                 }), more=FALSE, position=c(xmin=0.5, ymin=0, xmax=1, ymax=0.5))
  }
}

#' Compute the broad-sense heritability from a mixed model analysis with
#' ASReml.
#' 
#' The broad-sense heritability (\eqn{H^2}) is computed from the variance
#' components from a linear mixed model \describe{ 
#' \item{list(list("V_G"))}{the genetic variance components} 
#' \item{list(list("V_E"))}{the residual variance component} 
#' }
#' 
#' The two approaches for calculating \eqn{H^2}: \describe{ 
#' \item{which = }{\deqn{H^2 = \frac{V_G}{V_G + V_E}}}
#' \item{list("\"overall\"")}{ \deqn{H^2 = \frac{V_G}{V_G + V_E}}} 
#' \item{which = }{ \deqn{H^2 = \frac{V_G}{V_G + V_E/r}}}
#' \item{list("\"overall\"")}{ \deqn{H^2 = \frac{V_G}{V_G + V_E/r}}}
#' }
#' @param model An asreml-object, see \code{\link[asreml]{asreml}}.
#' @param which A character string/vector specifying which heritability measure
#' to compute, choose from \code{"line-mean"} and \code{"overall"}.  By default
#' both measures are computed.
#' @param genetic.factor (Optional) Character string specifying the name of the
#' genetic factor in the model. Default = \code{"GERMPLASM"}.
#' @param by.location Logical indicating whether to compute the heritability on
#' location-basis or across locations. Default = \code{FALSE}. 
#' @param data (Optional) The data frame used in the preceeding
#' \code{asreml}-call.
#' @param checks (Optional) The names of the check-lines in the
#' \code{genetic.factor}-column in \code{data}.
#' @return A data frame with the estimated genetic-, GxE- and residual standard
#' deviations and the computed heritability.
#' @author Ruud Derijcker, Katrien Baert
#' @export 
#' 
calculateH2 <- function(model, which = c("line-mean", "overall"), 
                        genetic.factor = "GERMPLASM",
                        by.location = FALSE, data = NULL, checks = NULL) {
  output <- list()
  varcomp <- summary(model)$varcomp
  if(by.location) {
    message("Location-based heritability is estimated only when gxe factor is included.")
    by.location <- FALSE
  }
  
  if(by.location) {
    locations <- unlist(strsplit(rownames(varcomp[grep("!variance", rownames(varcomp)), ]), split="_"))
    locations <- locations[grep("!", locations)]
    locations <- unlist(strsplit(locations, split="!"))
    locations <- locations[-grep("variance", locations)]
    if(length(locations) == 1 & locations[1] == "R") {
      by.location <- FALSE
    }
  }
  
  # Extract GENETIC variance
  var.gen <- varcomp[grep(paste(genetic.factor, "!", genetic.factor, ".var", sep = ""), 
                          rownames(varcomp)), "component"]
  if(length(var.gen) == 0) {
    # Try component modelled with giv-matrix
    var.gen <- varcomp[grep(paste("giv\\(", genetic.factor, "\\)\\.giv", sep = ""), 
                            rownames(varcomp)), "component"]
    
    if(length(var.gen) == 0)
      stop("The factor expressing genetic variation was not found.")
  }
  
  # Extract GxE variance
  var.gxe <- 0
  
  # Extract RESIDUAL variance (average if more than one)
  var.res <- varcomp[grep("!variance", rownames(varcomp)), "component"]
  if(length(var.res) > 0 & !by.location){
    var.res <- mean(var.res)
  }else if(length(var.res) != length(locations)){
    var.res <- mean(var.res)
  }
  
  # Prepare dataset
  if("line-mean" %in% which){
    if(is.null(data) || nrow(data) != length(model$residuals))
      stop("For computing the 'line-mean' heritability, please provide a the dataset you used for running asreml in the 'data' argument.")
    
    if(!is.null(checks)) {
      data[, genetic.factor] <- as.character(data[, genetic.factor])
      checks.in <- which(checks %in% data[, genetic.factor])
      if(length(checks.in) == 0){
        message("NOTE: None of the specified 'checks' is found in the ", genetic.factor, " column.")
        checks <- NULL
      } else if(length(checks.in) < length(checks)) {
        message("NOTE: Not all specified 'checks' are found in the ", genetic.factor, " column.")
        checks <- checks[checks.in]
      }
      data <- data[!data[, genetic.factor] %in% checks, ]
    }
  }  
  
  ## ACTUAL CALCULATIONS
  if("overall" %in% which) {
    hb <- var.gen/(var.gen + var.gxe + var.res)
    output[["overall"]] <- data.frame(STD.Gen = sqrt(var.gen), 
                                      STD.GxE = ifelse(var.gxe > 0, sqrt(var.gxe), NA), 
                                      R.GxE = NA,
                                      STD.Res = sqrt(var.res), 
                                      R.Res = NA,
                                      H2 = round(hb, 2))
    if(by.location) {
      oldColumns <- colnames(output[["overall"]])
      output[["overall"]]$LOCATION <- locations
      output[["overall"]] <- output[["overall"]][, c("LOCATION", oldColumns)]
    }
  }
  if("line-mean" %in% which) {
    data[, genetic.factor] <- factor(data[, genetic.factor])
    t.res <- table(table(data[, genetic.factor]))
    t.res <- t.res[names(t.res) != "0"]
    if(length(t.res) > 1) {
      n.res <- as.numeric(names(t.res[order(t.res, decreasing = TRUE)[1]]))
      message("NOTE: The design was not completely balanced, number of overall replicates choosen to be: ", n.res)
    } else {
      n.res <- as.numeric(names(t.res))
    }      
    
    n.gxe <- 1
    if(all(var.gxe != 0)) {
      gxe.factor1 <- gsub(genetic.factor, "", gxe.factor)
      gxe.factor1 <- gsub(":", "", gxe.factor1)
      data[, gxe.factor1] <- factor(data[, gxe.factor1])
      
      #			t.gxe <- table(table(data[, genetic.factor], data[, gxe.factor1]))
      t.gxe <- table(apply(table(data[, genetic.factor], data[, gxe.factor1]), 1, function(x) sum(x != 0)))
      t.gxe <- t.gxe[names(t.gxe) != "0"]
      if(length(t.res) > 1) {
        n.gxe <- as.numeric(names(t.gxe[order(t.gxe, decreasing = TRUE)[1]]))
        message("NOTE: The design was not completely balanced, number of GxE replicates choosen to be: ", n.gxe)
      } else {
        n.gxe <- as.numeric(names(t.gxe))
      }
    }
    if(!by.location) {
      hb <- var.gen/(var.gen + var.gxe/n.gxe + var.res/n.res)
      output[["line-mean"]] <- data.frame(STD.Gen = sqrt(var.gen), 
                                          STD.GxE = ifelse(var.gxe > 0, sqrt(var.gxe), NA), 
                                          R.GxE = ifelse(var.gxe > 0, n.gxe, NA),
                                          STD.Res = sqrt(var.res), 
                                          R.Res = n.res,
                                          H2 = round(hb, 2))
    } else {
      repByLocation <- function(x) {
        t <- table(x)
        t <- t[names(t) != 0]
        t1 <- as.numeric(names(t[order(t, decreasing=TRUE)][1]))
      }
      n.gxe <- 1
      n.res <- apply(table(data[, genetic.factor], data[, gxe.factor1]), 2, repByLocation)
      hb <- var.gen/(var.gen + var.gxe/n.gxe + var.res/n.res)
      output[["line-mean"]] <- data.frame(STD.Gen = sqrt(var.gen), 
                                          STD.GxE = ifelse(var.gxe > 0, sqrt(var.gxe), NA), 
                                          R.GxE = ifelse(var.gxe > 0, n.gxe, NA),
                                          STD.Res = sqrt(var.res), 
                                          R.Res = n.res,
                                          H2 = round(hb, 2))
      oldColumns <- colnames(output[["line-mean"]])
      output[["line-mean"]]$LOCATION <- locations
      output[["line-mean"]] <- output[["line-mean"]][, c("LOCATION", oldColumns)]
    }
    
  }
  return(list.to.dataframe(output, slotname = "Method")[, -1])
}