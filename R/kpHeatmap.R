#' kpHeatmap
#' 
#' @description 
#' 
#' Plots the given data as a heatmap along the genome
#' 
#' @details 
#'  
#' Given regions of the genome with a start, end and a value, draws a heatmap-like
#' representation, with the color of the region determined by its value. It is important to 
#' note that \code{kpHeatmap} will not extend the regions in any way, so if regions are not 
#' contiguous, they will appear as a series of rectangles and not as a continuous plot.
#' 
#' @usage kpHeatmap(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, autotrack=NULL, data.panel=1, colors=c("blue", "white", "yellow"), clipping=TRUE, ...)
#'  
#' @inheritParams kpPoints
#' @param x0 (numeric) the position (in base pairs) where the data region starts
#' @param x1 (numeric) the position (in base pairs) where the data region ends
#' @param colors    (colors) A set of color used to determine the color associated with each value. Internally, it uses \code{\link[grDevices]{colorRamp}}. (defaults to c("blue", "white", "yellow"))
#'     
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpLines}}
#' 
#' @examples
#'   
#' dd <- toGRanges(data.frame(chr="chr1", start=4980000*(0:49), end=4980000*(1:50)))
#' y <- sin(x=c(1:length(dd))/2)
#' 
#' kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#' 
#' kpLines(kp, dd, y=y, r0=0.4, r1=0.6, ymin=-1, ymax=1)
#' kpAxis(kp, r0=0.4, r1=0.6, ymin=-1, ymax=1, cex=0.5)
#' 
#' kpHeatmap(kp, dd, y=y, colors = c("red", "black", "green"), r0=0, r1=0.2)
#' kpHeatmap(kp, dd, y=y, colors = c("green", "black", "red"), r0=0.2, r1=0.4)
#' 
#' #or we can provide all data into a single GRanges object
#' mcols(dd) <- data.frame(y=y)
#' 
#' kpHeatmap(kp, dd, r0=0.6, r1=0.8)
#' #non-contiguous regions appear as solitary rectangles
#' kpHeatmap(kp, sample(x = dd, 10), r0=0.8, r1=1, color=c("orange", "black", "purple", "green"))
#' 
#' 
#'@export kpHeatmap

kpHeatmap <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, 
                      ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, 
                      autotrack=NULL, data.panel=1, 
                      colors=c("blue", "white", "yellow"), clipping=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
 
  #Manually process "y" since in this track it is not affected by r0, r1, autotrack, etc...
  if(!is.null(data)) {
    if(is.null(y)) {
      if("value" %in% names(mcols(data))) {
        y <- data$value
      } else {
        if("y" %in% names(mcols(data))) {
          y <- data$y
        } else {
          stop("No y value specified. It is needed to provide ymax or a column named 'y' in data")
        }
      }
    }
  } 
  
  if(is.null(ymin)) ymin <- min(y)
  if(is.null(ymax)) ymax <- max(y)
  
  if(ymin == ymax) {
    ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
    ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  }
      
  #Standardize the values to the [0,1] range to plot it with colorRamp
  y <- y - ymin
  y <- y/(ymax - ymin)
  
  
  #Create the colorRamp
  cr <- grDevices::colorRamp(colors=colors)
  
  invisible(kpRect(karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=0, y1=1,
                   ymin=0, ymax=1, r0=r0, r1=r1, autotrack=autotrack,
                   col=grDevices::rgb(cr(y), max=255), border=NA,
                   data.panel=data.panel, clipping=clipping, ...))
  
  
 
}
