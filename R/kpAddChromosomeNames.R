#' kpAddChromosomeNames
#' 
#' @description 
#' 
#' Plots the chromosome names in the karyoplot
#' 
#' @details 
#' 
#' Given a KaryoPlot object, plot the names of the depicted chromosomes. This 
#' function is usually automatically called by plotKaryotype unless 
#' \code{labels.plotter} is NULL.
#' 
#' @usage kpAddChromosomeNames(karyoplot, names=NULL, xoffset=0, yoffset=0, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param names      (character vector) the names to use for the chromosomes. If NULL, the chromosome names in the original genome will be used. (defaults to NULL)
#' @param xoffset    (numeric) a number of units to move the the chromosome names on the x axis with respect to their standard position (defaults to 0)
#' @param yoffset    (numeric) a number of units to move the the chromosome names on the y axis with respect to their standard position (defaults to 0)
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{getChromosomeNamesBoundingBox}}
#' 
#' @examples
#'
#' kp <- plotKaryotype(labels.plotter = NULL)
#' kpAddChromosomeNames(kp, col="red", srt=30)
#'  
#' @export kpAddChromosomeNames
#' 

kpAddChromosomeNames <- function(karyoplot, chr.names=NULL, xoffset=0, yoffset=0, ...) {
  #Validate parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(is.null(chr.names)) chr.names <- karyoplot$chromosomes
  if(length(chr.names)==0) stop("In kpAddChromosomeNames: chr.names must have at least one element.")
  if(!all(methods::is(chr.names, "character"))) stop("In kpAddChromosomeNames: all elements of chr.names must be characters.")

  #Begin plotting
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  bb <- getChromosomeNamesBoundingBox(karyoplot)
  
  x <- (bb$x0+bb$x1)/2 + xoffset
  y <- (bb$y0+bb$y1)/2 + yoffset
  
  graphics::text(x=x, y=y, labels=chr.names, ...)
  
  invisible(karyoplot)
}

