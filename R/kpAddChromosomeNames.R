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
#' @usage kpAddChromosomeNames(karyoplot, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
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

kpAddChromosomeNames <- function(karyoplot, ...) {
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  bb <- getChromosomeNamesBoundingBox(karyoplot)
  
  chr.labels <- karyoplot$chromosomes

  x <- (bb$x0+bb$x1)/2
  y <- (bb$y0+bb$y1)/2
  
  graphics::text(x=x, y=y, labels=chr.labels, ...)
  
  invisible(karyoplot)
}

