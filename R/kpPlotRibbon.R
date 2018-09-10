#' kpPlotRibbon
#' 
#' @description 
#' 
#' A variable width ribbon
#' 
#' @details 
#'  
#' \code{kpPlotRibbon} plots a variable witdh ribbon along the genome. It can be used,
#' for example, to plot the sd region around a line representing a mean. It can also 
#' be used as a replacement for \code{\link{kpBars}} creating a smoother plot without 
#' the the actual individual bars. \code{kpPlotRibbon} has three additional parameters
#' controlling the smoothing of the lines and their colors.
#' 
#' @usage kpPlotRibbon(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col="gray80", border=NULL, clipping=TRUE, ...)
#'  
#' @inheritParams kpRect 
#' @param col  (color) The background color to plot. If NULL, it will be a lighter version of 'border' or 'black' if border is null. (Defaults to "gray80")
#' @param border (color) The color to use to plot the borders of the bars. If NULL, it will be a darker version of 'col'. If NA, no border will be plotted. (Defaults to NULL)
# @param smooth A boolean indicating if the ribbon should be smoothed (with a spline approximation) before plotting (defaults to FALSE)
#'     
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpBars}}, \code{\link{kpLines}}
#' 
#' @examples
#' 
#' 
#' set.seed(1000)
#' 
#' data <- toGRanges(data.frame(chr="chr1", start=1e6*(0:239), end=1e6*(1:240)))
#' y <- ((sin(start(data))/5 + rnorm(n=24, mean=0, sd=0.1))/5)+0.5
#'  
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes="chr1")
#' 
#' kpPlotRibbon(kp, data, y0=y-0.3, y1=y+0.3, border="red", col=lighter("red"))
#' kpPlotRibbon(kp, data, y0=y-0.1, y1=y+0.1, border="blue", col=lighter("blue"))
#' kpLines(kp, data, y=y, col="green")
#' kpPlotRibbon(kp, data, y0=0.5+(y-min(y)), y1=0.5-(y-min(y)), data.panel=2)
#' 
#' @export kpPlotRibbon


# smooth=FALSE,
kpPlotRibbon <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col="gray80", border=NULL, clipping=TRUE, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotRibbon: 'karyoplot' must be a valid 'KaryoPlot' object"))
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #If y0 is not specified with any of the valid methods, set it to the min of the data.panel
  if(is.null(y0)) {
    if(is.null(data)) {
      y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
    } else {
      if(!("y0" %in% names(mcols(data)))) {
        y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
      }
    }
  }
  
  #Specify the missing colors if possible
  if(is.null(border) & !is.null(col) & !is.na(col)) {
    border=darker(col, amount = 100)
  }
  if(is.null(col) & !is.null(border) & !is.na(border)) {
    col=lighter(border)
  }
  if(is.na(col) & is.null(border)) {
    border <- "black"
  }
  
  pp <- prepareParameters4("kpPlotRibbon", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1,
                           y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, 
                           autotrack=autotrack, data.panel=data.panel, ...)
  
  ccf <- karyoplot$coord.change.function
  
  x0plot <- ccf(chr=pp$chr, x=pp$x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=pp$chr, x=pp$x1, data.panel=data.panel)$x
  #TODO: Add a parameter to specify center (default), start (xplot=x0plot) or end(xplot=x1plot)
  xplot <- (x0plot+x1plot)/2
  y0plot <- ccf(chr=pp$chr, y=pp$y0, data.panel=data.panel)$y
  y1plot <- ccf(chr=pp$chr, y=pp$y1, data.panel=data.panel)$y
  
  for(chr.name in karyoplot$chromosomes) {
    
    in.chr <- pp$chr==chr.name
    chr.chr <- pp$chr[in.chr]
    x.chr <- xplot[in.chr]
    y0.chr <- y0plot[in.chr]
    y1.chr <- y1plot[in.chr]
    
    if(karyoplot$zoom==TRUE) {
      if(clipping==TRUE) {
        dpbb <- karyoplot$getDataPanelBoundingBox(data.panel)
        graphics::clip(x1 = dpbb$xleft, x2 = dpbb$xright, y1 = dpbb$ybottom, y2=dpbb$ytop)
      }
    }
    
    
    
    graphics::polygon(x=c(x.chr, rev(x.chr)), y=c(y0.chr, rev(y1.chr)), col=col, border=NA, ...)
    # if(smooth==TRUE) {
    #   graphics::lines(smooth.spline(x = x.chr, y=y0.chr, df = length(x.chr)), col=border, ...)
    #   graphics::lines(smooth.spline(x = x.chr, y=y1.chr, df = length(x.chr)), col=border, ...)
    # } else {
      graphics::lines(x=x.chr, y=y0.chr, col=border, ...)
      graphics::lines(x=x.chr, y=y1.chr, col=border, ...)
    # }
  }
  
  invisible(karyoplot)
}

