#' getMainTitleBoundingBox
#' 
#' @description 
#' 
#' Return the regions where the chromosome names should be placed
#' 
#' @details 
#' 
#' Given a KaryoPlot object, return the regions where the main plot
#' should be placed. The position will depend on the plot type used.
#' 
#' @note In general, this function is automatically called by karyoploteR
#' and the user never needs to call it. 
#' 
#' @usage getMainTitleBoundingBox(karyoplot)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' 
#' @return
#' Returns a list with four elements (x0, x1, y0 and y1), each of them an 
#' integer with the coordinates for the main title
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpAddMainTitle}}
#' 
#' @examples
#'
#' kp <- plotKaryotype()
#' bb <- getMainTitleBoundingBox(kp)
#'  
#' @export getMainTitleBoundingBox
#' 


getMainTitleBoundingBox <- function(karyoplot) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(karyoplot$plot.type %in% c(1,2,4,6,7)) {
    y0 <- karyoplot$plot$ymax - karyoplot$plot.params$topmargin
    y1 <- karyoplot$plot$ymax
    x0 <- karyoplot$plot$xmin
    x1 <- karyoplot$plot$xmax
  
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
  if(karyoplot$plot.type %in% c(3, 5)) {
    y0 <- karyoplot$plot$ymax - karyoplot$plot.params$topmargin/2 #The bottom part of the top margin is used by the chromosome names
    y1 <- karyoplot$plot$ymax
    x0 <- karyoplot$plot$xmin
    x1 <- karyoplot$plot$xmax
    return(list(x0=x0, x1=x1, y0=y0, y1=y1))  
  }
 
  
  #Else
  stop("Unknown plot type: ", karyoplot$plot.type)
}
