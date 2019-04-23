#' kpPlotNames
#' 
#' @description 
#' 
#' Plots text labels with positioning relative to rectangles along the genome. 
#' 
#' @details 
#'  
#'  This is a simple wrapper around \code{\link{kpText}} that positions the 
#'  text relative to the rectangles defined by its arguments. They may be
#'  used to name or label different graphical elements in the plot.
#'  The rectangles may be specified as in \code{\link{kpRect}} the relative 
#'  positions accepted are: "left", "right", "top", "bottom", "center". 
#'  It is possible to specify and empty label (\code{labels=""}) to leave an
#'  element without name.
#'
#' @usage kpPlotNames(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=NULL, y1=NULL, labels=NULL, position="left", ymax=NULL, ymin=NULL, r0=NULL, r1=NULL,  data.panel=1, clipping=TRUE, ...)
#'
#' @inheritParams kpRect
#' @param position (character) The position of the text relative to the rectangle. Can be "left", "right", "top", "bottom" or "center". Defaults to "left".
#' @param labels  (character) The labels to use in the plot. They will be associated to the rectangles by its order and recycled as needed. 
#'  
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{kpText}}, \code{\link{kpRect}}
#' 
#' @examples
#'  
#'  
#'  regs <- toGRanges(data.frame(chr=c("chr1", "chr1", "chr1"),
#'                   start=c(20e6, 100e6, 200e6),
#'                   end=c(40e6, 170e6, 210e6),
#'                   y0=c(0.1, 0.5, 0.7),
#'                   y1=c(0.5, 0.6, 0.95)))
#'  
#'  kp <- plotKaryotype(genome="hg19", chromosomes="chr1")
#'  kpRect(kp, data=regs)  
#'  
#'  kpPlotNames(kp, data=regs, labels=c("R1", "R2", "R3"))
#'  kpPlotNames(kp, data=regs, labels=c("R1", "R2", "R3"), position="top", cex=2)    
#'  kpPlotNames(kp, data=regs, labels=c("R1", "", "R3"), position="right", col="red")
#'  kpPlotNames(kp, data=regs, labels="bottom", position="bottom", col=rainbow(3))
#'  kpPlotNames(kp, data=regs, labels="o", position="center", col=rainbow(3), cex=1)                
#'  
#'@export kpPlotNames


kpPlotNames <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=NULL, y1=NULL, 
                        labels=NULL, position="left",
                        ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, 
                        data.panel=1, clipping=TRUE, ...) {

  
  #karyoplot
    if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #position
  position <- match.arg(position, c("left", "right", "top", "bottom", "center"))
  
  if(length(data)==0 && length(chr)==0) return(invisible(karyoplot))

  #Note: we use r0=0 and r1=1 (and ymin=0 and ymax=1) so we only use data input normalization (data, chr, etc...) but not normalization of r0, ymins, etc...  
  pp <- prepareParameters4("kpPlotNames", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1,
                           y0=y0, y1=y1, ymin=0, ymax=1, r0=0, r1=1, 
                           data.panel=data.panel, ...)
  
  #if there's nothing to plot, return
  if(length(pp$chr)==0) {
    invisible(karyoplot)
  }
  
  #Now decide how to plot (with respect to the rectangles), and call kpText with the appropiate parameters
  switch(position,  
    left=kpText(karyoplot, chr=pp$chr, x=pp$x0, y=pp$y0+(pp$y1-pp$y0)/2, labels=labels, pos=2, ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...),
    right=kpText(karyoplot, chr=pp$chr, x=pp$x1, y=pp$y0+(pp$y1-pp$y0)/2, labels=labels, pos=4, ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...),
    top=kpText(karyoplot, chr=pp$chr, x=pp$x0+(pp$x1-pp$x0)/2, y=pp$y1, labels=labels, pos=3,  ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...),
    bottom=kpText(karyoplot, chr=pp$chr, x=pp$x0+(pp$x1-pp$x0)/2, y=pp$y0, labels=labels, pos=1,  ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...),
    center=kpText(karyoplot, chr=pp$chr, x=pp$x0+(pp$x1-pp$x0)/2, y=pp$y0+(pp$y1-pp$y0)/2, labels=labels, ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...)
  )  

  invisible(karyoplot)
}
