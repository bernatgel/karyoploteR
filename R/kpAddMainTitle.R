#' kpAddMainTitle
#' 
#' @description 
#' 
#' Plots the chromosome names in the karyoplot
#' 
#' @details 
#' 
#' Given a KaryoPlot object and a character string, plot the character strings
#' as the main title of the plot. This function is usually automatically
#' called by plotKaryotype unless.
#' 
#' @usage kpAddMainTitle(karyoplot, main=NULL, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param main (character) the main title of the plot
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{getMainTitleBoundingBox}}
#' 
#' @examples
#'
#' kp <- plotKaryotype(labels.plotter = NULL)
#' kpAddMainTitle(kp, col="red", srt=30, cex=0.8)
#'  
#' @export kpAddMainTitle
#' 

kpAddMainTitle <- function(karyoplot, main=NULL, ...) {
  
  if(!is.null(main)) {
    karyoplot$beginKpPlot()
    on.exit(karyoplot$endKpPlot())
    
    bb <- getMainTitleBoundingBox(karyoplot)
    
    x <- (bb$x0+bb$x1)/2
    y <- (bb$y0+bb$y1)/2
    
    graphics::text(x=x, y=y, labels=main, ...)
    
  }
  
  invisible(karyoplot)
}

