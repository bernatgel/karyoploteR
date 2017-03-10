#' getChromosomeNamesBoundingBox
#' 
#' @description 
#' 
#' Return the regions where the chromosome names should be placed
#' 
#' @details 
#' 
#' Given a KaryoPlot object, return the regions where the chromosome labels 
#' should be placed. The positions will depend on the plot type used.
#' 
#' @note In general, this function is automatically called by karyoploteR
#' and the user never needs to call it. 
#' 
#' @usage getChromosomeNamesBoundingBox(karyoplot)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' 
#' @return
#' Returns a list with four elements (x0, x1, y0 and y1), each of them a names
#' vector of integers with one coordinatefor every chromosome in the plot.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{plotChromosomeNames}}
#' 
#' @examples
#'
#' kp <- plotKaryotype(ideogram.plotter = NULL)
#' bb <- getChromosomeNamesBoundingBox(kp)
#'  
#' @export getChromosomeNamesBoundingBox
#' 


getChromosomeNamesBoundingBox <- function(karyoplot) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(kp$plot.type %in% c(1,2)) {
    chr.labels <- karyoplot$chromosomes
    
    y0 <- karyoplot$ideogram.mid(chr=chr.labels) - kp$plot.params$ideogramheight
    y1 <- karyoplot$ideogram.mid(chr=chr.labels) + kp$plot.params$ideogramheight
    x0 <- setNames(rep(0, length(chr.labels)), chr.labels)
    x1 <- setNames(rep(karyoplot$plot.params$leftmargin, length(chr.labels)), chr.labels)
  
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
  
  #Else
  stop("Unknown plot type")
}
