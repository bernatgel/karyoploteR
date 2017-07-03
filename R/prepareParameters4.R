#' prepareParameters4
#' 
#' @description 
#' Prepare and normalize the parameters for functions with x0, x1 and y0, y1 parameters
#'  
#' @details 
#' This function prepares and normalizes the parameters for plotting functions  
#' with x0, x1, y0 and y1 parameters (as opposed to x and y) so functions can
#' offer a richer interface while internally dealing only with standard and 
#' simple code. It extracts the 
#' positions from \code{data} if available and applies the \code{r0} and 
#' \code{r1} scaling. It returns the ready to plot values in a list with
#' only \code{chr}, \code{x0}, \code{x1}, \code{y0} and \code{y1}. 
#' Individual parameters (\code{chr}, \code{x0}, \code{x1}, \code{y0} and 
#' All parameters are interpreted and used as explained in \code{\link{kpRect}}.
#' It also filters out any data points corresponding to chromosomes not present
#' in the current karyoplot.
#'  
#' @note This function is only useful when creating custom plotting functions. 
#' It is not intended to the general user.
#' 
#' @note For detailed documentation on the parameters, see \code{\link{kpRect}}
#'  
#' @usage prepareParameters4(function.name, karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, filter.data=TRUE, ...)
#'  
#' @param function.name (character) The name of the function calling \code{prepareParameters4}. Only user for error reporting.
#' @param karyoplot (KaryoPlot) A karyoplot object.
#' @param data A GRanges
#' @param chr A character representing the chromosome names.
#' @param x0 The position in the chromosome in number of bases.
#' @param x1 The position in the chromosome in number of bases.
#' @param y0 The value to be plotted.
#' @param y1 The value to be plotted.
#' @param ymax The maximum value of y
#' @param ymin The minimum value of y
#' @param r0 The start of the range to use for plotting
#' @param r1 The end of the range to use for plotting
#' @param data.panel The data panel to use
#' @param filter.data A boolean indicating if data should be filtered so only data in visible chromosomes is kept. (defaults to TRUE, filter data)
#' @param ... Any additional parameter
#'
#' @return 
#' A list with five values: \code{chr}, \code{x0}, \code{x1}, \code{y0} and \code{y1}. Each of them 
#' a vector of the same length with the normalized values to plot.
#'
#'
#' @seealso \code{\link{kpRect}}
#' 
#' @examples
#' 
#' kp <- plotKaryotype()
#' prepareParameters4("TestFunc", kp, data=NULL, chr="chr1", x0=c(10, 20, 30), x1=c(20, 30, 40), y0=c(0, 1, 2), y1=c(0.5, 1.5, 3), r0=0, r1=0.5, ymin=0, ymax=3)
#' 
#'  
#' @export prepareParameters4
#' 


prepareParameters4 <- function(function.name, karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, filter.data=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In ", function.name, ": 'karyoplot' must be a valid 'KaryoPlot' object"))
  
  #if null or NA, get the r0 and r1 and ymin-ymax from the plot params
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.na(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.na(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.na(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.na(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  
  
  if(!is.null(data)) {
    if(all(!is.na(data))) {
      if(is.null(chr) || all(is.na(chr))) {
        chr <- as.character(seqnames(data))
      }
      if(is.null(x0) || all(is.na(x0))) {
        x0 <- start(data)
      }
      if(is.null(x1) || all(is.na(x1))) {
        x1 <- end(data)
      }
      
      if(is.null(y0) || all(is.na(y0))) {
        if("y0" %in% names(mcols(data))) {
          y0 <- data$y0
        } else {
          stop("No y0 value specified. Parameter y0 or a column named 'y0' in data must be provided")
        }
      }
      if(is.null(y1) || all(is.na(y1))) {
        if("y1" %in% names(mcols(data))) {
          y1 <- data$y1
        } else {
          stop("No y1 value specified. Parameter y1 or a column named 'y1' in data must be provided")
        }
      }
    }
  } 
    
  if(is.null(chr)) stop("chr must be specified, either by the 'chr' parameter or by providing a 'data' object")
  #if(any(is.na(chr))) stop("chr cannot be NA")   
  if(is.null(x0)) stop("x0 must be specified, either by the 'x0' parameter or by providing a 'data' object")
  #if(any(is.na(x0))) stop("x0 cannot be NA")   
  if(is.null(y0)) stop("y0 must be specified, either by the 'y0' parameter or by providing a 'data' object with a column named 'y0'")
  #if(any(is.na(y0))) stop("y0 cannot be NA")   
  if(is.null(x1)) stop("x1 must be specified, either by the 'x1' parameter or by providing a 'data' object")
  #if(any(is.na(x1))) stop("x1 cannot be NA")   
  if(is.null(y1)) stop("y1 must be specified, either by the 'y1' parameter or by providing a 'data' object with a column named 'y1'")
  #if(any(is.na(y1))) stop("y1 cannot be NA")   
  
  #transform chr to a character
  chr <- as.character(chr)

  #Check the classes
  if(!base::is.character(chr)) stop("chr must be a character vector")   
  if(!base::is.numeric(x0)) stop("x0 must be a numeric vector")
  if(!base::is.numeric(x1)) stop("x1 must be a numeric vector")
  if(!base::is.numeric(y0)) stop("y0 must be a numeric vector")
  if(!base::is.numeric(y1)) stop("y1 must be a numeric vector")   
  
    
  #Recicle any values as needed
  chr <- recycle.first(chr, x0, x1, y0, y1)
  x0 <- recycle.first(x0, chr, x1, y0, y1)
  x1 <- recycle.first(x1, chr, x0, y0, y1)
  y0 <- recycle.first(y0, chr, x0, x1, y1)
  y1 <- recycle.first(y1, chr, x0, x1, y0)

  
  #Scale it with ymin and ymax
  y0 <- (y0 - ymin)/(ymax - ymin)
  y1 <- (y1 - ymin)/(ymax - ymin)
  
  #scale y to fit in the [r0, r1] range
  y0 <- (y0*(r1-r0))+r0
  y1 <- (y1*(r1-r0))+r0
  
  if(filter.data) {
    in.visible.chrs <- chr %in% karyoplot$chromosomes
    chr <- chr[in.visible.chrs]
    x0 <- x0[in.visible.chrs]
    x1 <- x1[in.visible.chrs]
    y0 <- y0[in.visible.chrs]
    y1 <- y1[in.visible.chrs]
  }
  
  
  return(list(chr=chr, x0=x0, x1=x1, y0=y0, y1=y1))
  
}
