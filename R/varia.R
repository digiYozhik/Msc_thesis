#' Convert data frame to list
#' 
#'@param x data frame
#'@param by character specifying column name where levels of factor are used to 
#'make the elements of a list
#'@return list with elements named after levels of the column specified by the 
#'by-argument.
#'@export
#'@author Ruud Derijcker
#'
dataframe.to.list <- function (x, by) {
  if (by == "rownames") 
    out <- split(x, rownames(x))
  else if (by %in% colnames(x)) 
    out <- split(x, x[, by])
  else stop("Invalid 'by'-argument.")
  return(out)
}

#' Transform a list of data frames to one data frame
#' 
#' @param x A \code{list} with a data frame with equal column names in each
#' slot and marker names as the row names.
#' @param slotname The column name for the original slot names, default =
#' \code{"Slot"}.
#' @param rowname The column name for the original row names, default =
#' \code{"Row"}.
#' @return A data frame.
#' @author Ruud Derijcker
#' @export
#' 
list.to.dataframe <- function(x, slotname = "Slot", rowname = "Row"){
  report <- NULL
  if(is.null(names(x)))
    names(x) <- 1:length(x)
  for(i in names(x)){
    temp <- data.frame(x[[i]])
    temp.colnames <- colnames(temp)
    temp$id <- i
    temp$row <- rownames(temp)
    rownames(temp) <- NULL
    report <- rbind(report, temp)
  }
  names(report)[names(report) == "id"] <- slotname
  names(report)[names(report) == "row"] <- rowname
  return(report[, c(rowname, slotname, temp.colnames)])
}