#' kpAbline
#' 
#' @description 
#' 
#' This is the \code{KaryoploteR} version of the \code{\link[graphics]{abline}} function to add horizontal or vertical lines to the plot.
#' 
#' 
#' @details 
#'  
#' As with all other "base-inspired" low-level plotting functions in karyoploteR, the function has been designed to accept mostly
#' the same parameters as the base one (see the package vignette for more information}. In this case, however, the interface has been reduced and it is only possible to plot 
#' vertical and horizontal lines and it's not possible to provide an intercept and slope. In addition, the function accepts
#' graphical parameters that are valid for the base function \code{\link[graphics]{segments}}.
#' 
#' 
#' @usage kpAbline(karyoplot, chr=NULL, h=NULL, v=NULL, ymin=NULL, ymax=NULL, data.panel=1,  r0=NULL, r1=NULL, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param chr    (a charecter vector) A vector of chromosome names specifying the chromosomes where the lines will be plotted. If NULL, the lines will be plotted in all chromosomes. (defaults to NULL)
#' @param h    (a numeric vector) A numeric vector with the heights where the horizontal lines will be plotted. If \code{h} is NULL, no horizontal lines will be plotted. (defaults to NULL)    
#' @param v    (a numeric vector) A numeric vector with the positions (in base pairs) where the vertical lines will be plotted. If \code{v} is NULL, no vertical lines will be plotted. (defaults to NULL)
#' @param ymin    (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL)
#' @param ymax    (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#'  
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpSegments}}, \code{\link{kpLines}}
#' 
#' @examples
#' 
#'   
#' 
#' @export kpAbline


kpAbline <- function(karyoplot, chr=NULL, h=NULL, v=NULL, ymin=NULL, ymax=NULL, data.panel=1,  r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
 
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  ccf <- karyoplot$coord.change.function
  
  if(is.null(chr)) { #if chr is not specified, plot the line in all chromosomes
    chr <- as.character(seqnames(karyoplot$genome)) 
  } 
    
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  
  if(!is.null(h)) {
    names(karyoplot$genome) <- as.character(seqnames(karyoplot$genome))
    chr <- rep(chr, each=length(h))
    x0 <- start(karyoplot$genome[chr])
    x1 <- end(karyoplot$genome[chr])
    
    kpSegments(karyoplot=karyoplot, chr=chr, x0=x0, x1=x1, y0=h, y1=h, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)  
  }    
   
  if(!is.null(v)) {
    y0 <- ymin 
    y1 <- ymax
    
    kpSegments(karyoplot=karyoplot, chr=chr, x0=v, x1=v, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)  
  }
  
  return(karyoplot)
        
}
