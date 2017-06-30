#' prepareParameters2
#' 
#' @description 
#' Prepare and normalize the parameters for functions with x and y parameters
#'  
#' @details 
#' This function prepares and normalizes the parameters for plotting functions  
#' with x and y parameters (as opposed to x0, x1, y0 and y1) so functions can
#' offer a richer interface while internally dealing only with standard and 
#' simple code. It extracts the 
#' positions from \code{data} if available and applies the \code{r0} and 
#' \code{r1} scaling. It returns the ready to plot values in a list with
#' only \code{chr}, \code{x} and \code{y}. All parameters are interpreted and
#'  used as explained in \code{\link{kpPoints}}. It also filters out any data 
#'  points corresponding to chromosomes not present in the current karyoplot.
#'  
#' @note This function is only useful when creating custom plotting functions. 
#' It is not intended to the general user.
#' 
#' @note For detailed documentation on the parameters, see \code{\link{kpPoints}}
#'  
#' @usage prepareParameters2(function.name, karyoplot, data, chr, x, y, ymax, ymin, r0, r1, data.panel, filter.data=TRUE, ...)
#'  
#' @param function.name (character) The name of the function calling \code{prepareParameters2}. Only user for error reporting.
#' @param karyoplot (KaryoPlot) A karyoplot object.
#' @param data A GRanges
#' @param chr A character representing the chromosome names.
#' @param x The position in the chromosome in number of bases.
#' @param y The value to be plotted.
#' @param ymax The maximum value of y
#' @param ymin The minimum value of y
#' @param r0 The start of the range to use for plotting
#' @param r1 The end of the range to use for plotting
#' @param data.panel The data panel to use
#' @param filter.data A boolean indicating if data should be filtered so only data in visible chromosomes is kept. (defaults to TRUE, filter data)
#' @param ... Any additional parameter
#'
#' @return 
#' A list with three values: \code{chr}, \code{x} and \code{y}. Each of them 
#' a vector of the same length with the normalized values to plot.
#'
#'
#' @seealso \code{\link{kpPoints}}
#' 
#' @examples
#' 
#' kp <- plotKaryotype()
#' prepareParameters2("TestFunc", kp, data=NULL, chr="chr1", x=c(10, 20, 30), y=c(0, 1, 2), r0=0, r1=0.5, ymin=0, ymax=2)
#' 
#'  
#' @export prepareParameters2
#' 


prepareParameters2 <- function(function.name, karyoplot, data, chr, x, y, ymax, ymin, r0, r1, data.panel, filter.data=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In ", function.name, ": 'karyoplot' must be a valid 'KaryoPlot' object"))
    
  #if null or NA, get the r0 and r1 and ymin-ymax from the plot params
  if(is.null(r0) | is.na(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1) | is.na(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(ymin) | is.na(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax) | is.na(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x <- start(data) + (end(data) - start(data))/2 #Get the midpoints of the regions
    
    if(is.null(y)) {
      if("y" %in% names(mcols(data))) {
        y <- data$y
      } else {
        if("value" %in% names(mcols(data))) {
          y <- data$value
        } else {
          stop("No y value specified. Parameter y or a column named 'y' or 'value' in data must be provided")  
        }
      }
    } 
  } 
  
  if(is.null(chr) | is.na(chr)) stop("chr must be specified, either by the 'chr' parameter or by providing a 'data' object")
  if(is.null(x) | is.na(x)) stop("x must be specified, either by the 'x' parameter or by providing a 'data' object")
  if(is.null(y) | is.na(y)) stop("y must be specified, either by the 'y' parameter or by providing a 'data' object with a column names 'value' or 'y'")
  
  #transform chr to a character
  chr <- as.character(chr)
  
  #Recicle any values as needed
  chr <- recycle.first(chr, x, y)
  x <- recycle.first(x, chr, y)
  y <- recycle.first(y, chr, x)
  
  #Scale it with ymin and ymax
  y <- (y - ymin)/(ymax - ymin)
  
  
  #scale y to fit in the [r0, r1] range
  y <- (y*(r1-r0))+r0
  
  if(filter.data) {
    in.visible.chrs <- chr %in% karyoplot$chromosomes
    chr <- chr[in.visible.chrs]
    x <- x[in.visible.chrs]
    y <- y[in.visible.chrs]
  }
    
  return(list(chr=chr, x=x, y=y))
  
}
