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
#' Returns a list with four elements (x0, x1, y0 and y1), each of them a named
#' vector of integers with one coordinatefor every chromosome in the plot.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpAddChromosomeNames}}
#' 
#' @examples
#'
#' kp <- plotKaryotype()
#' bb <- getChromosomeNamesBoundingBox(kp)
#'  
#' @export getChromosomeNamesBoundingBox
#' 


getChromosomeNamesBoundingBox <- function(karyoplot) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(karyoplot$plot.type %in% c(1,2)) {
    chr.labels <- karyoplot$chromosomes
    
    y0 <- karyoplot$ideogram.mid(chr=chr.labels) - karyoplot$plot.params$ideogramheight
    y1 <- karyoplot$ideogram.mid(chr=chr.labels) + karyoplot$plot.params$ideogramheight
    x0 <- stats::setNames(rep(0, length(chr.labels)), chr.labels)
    x1 <- stats::setNames(rep(karyoplot$plot.params$leftmargin, length(chr.labels)), chr.labels)
  
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
  if(karyoplot$plot.type %in% c(3, 5)) {
    #position the labels centered in the chromosome and in the top margin
    chr.labels <- karyoplot$chromosomes
    
    chr.lens <- stats::setNames(as.numeric(end(karyoplot$plot.region) - start(karyoplot$plot.region)), chr.labels)
    
    y1 <- stats::setNames(rep(karyoplot$plot$ymax - karyoplot$plot.params$topmargin/2, length(chr.labels)), chr.labels) #the top half of the top margin is reserved for the main title
    y0 <- stats::setNames(rep(karyoplot$plot$ymax - karyoplot$plot.params$topmargin, length(chr.labels)), chr.labels)
    
    #use the coordinates change function to get the x positioning
    ccf <- karyoplot$coord.change.function
    x0 <- ccf(chr=chr.labels, x=start(karyoplot$plot.region), data.panel=1)$x
    x1 <- ccf(chr=chr.labels, x=end(karyoplot$plot.region), data.panel=1)$x
    
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
  if(karyoplot$plot.type %in% c(4)) {
    #position the labels centered in the chromosome and in the bottom margin
    chr.labels <- karyoplot$chromosomes
    
    chr.lens <- stats::setNames(as.numeric(end(karyoplot$plot.region) - start(karyoplot$plot.region)), chr.labels)
    
    y1 <- stats::setNames(rep(karyoplot$plot.params$bottommargin, length(chr.labels)), chr.labels)
    y0 <- stats::setNames(rep(0, length(chr.labels)), chr.labels)
    
    #use the coordinates change function to get the x positioning
    ccf <- karyoplot$coord.change.function
    x0 <- ccf(chr=chr.labels, x=start(karyoplot$plot.region), data.panel=1)$x
    x1 <- ccf(chr=chr.labels, x=end(karyoplot$plot.region), data.panel=1)$x
    
     
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
  
  #Else
  stop("Unknown plot type")
}
