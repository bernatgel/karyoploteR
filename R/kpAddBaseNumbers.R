#' kpAddBaseNumbers
#' 
#' @description 
#' 
#' Plots the base numbers along the chromosome ideograms
#' 
#' @details 
#'  
#' This function can be used to add the base numbers scale to the chromosome ideograms.
#' The base numbers and ticks witll be drawn next to the ideograms and not on a separate
#' independent x axis. It is possible to control the number and position of the tick
#' marks and labels
#' 
#' @usage kpAddBaseNumbers(karyoplot, tick.dist=20000000, tick.len=5, minor.ticks=TRUE, minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, tick.col=NULL, minor.tick.col=NULL, ...)
#' 
#' @param karyoplot  (karyoplot object) A valid karyoplot object created by a call to \code{\link{plotKaryotype}}
#' @param tick.dist  (numeric) The distance between the major numbered tick marks in bases
#' @param tick.len  (numeric) The length of the major tick marks in plot coordinates
#' @param minor.ticks (boolean) Whether to add unlabeled minor ticks between the major ticks
#' @param minor.tick.dist   (numeric) The distance between the minor ticks in bases
#' @param minor.tick.len   (numeric) The length of the minor tick marks in plot coordinates
#' @param cex  (numeric) The cex parameter for the major ticks label
#' @param tick.col   (color) If specified, the color to plot the major ticks. Otherwise the default color or, if given, the col parameter will be used. (Defaults to NULL)
#' @param minor.tick.col  (color) If specified, the color to plot the minor ticks. Otherwise the default color or, if given, the col parameter will be used. (Defaults to NULL)
#' @param ...  Any other parameter to be passed to internal function calls. Specially useful for graphic parameters.
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}
#' 
#' @examples
#'  
#' 
#' kp <- plotKaryotype()
#' kpAddBaseNumbers(kp)
#' 
#' kp <- plotKaryotype(chromosomes="chr17")
#' kpAddBaseNumbers(kp, tick.dist=10000000, minor.tick.dist=1000000)
#' 
#' 
#'  
#' @export kpAddBaseNumbers
#' 


kpAddBaseNumbers <- function(karyoplot, tick.dist=20000000, tick.len=5, minor.ticks=TRUE, 
                            minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, 
                            tick.col=NULL, minor.tick.col=NULL, ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  
  toLabel <- function(n) {
    if(abs(n) < 1000) return(as.character(n))
    if(abs(n) < 1000000) return(paste0(as.character(round(n/10)/100), "")) #Kb
    return(paste0(as.character(round(n/10000)/100), "")) #Mb
  }
  
  old.scipen <- options("scipen")
  options(scipen=999)
  on.exit(options(scipen=old.scipen), add=TRUE)
  
  #For every chromsome
  for(chr.name in as.character(seqnames(karyoplot$genome))) {

    chr <- karyoplot$genome[as.character(seqnames(karyoplot$genome)) == chr.name]
    #Major ticks
      num.ticks <- width(chr)/tick.dist + 1
      tick.pos <- start(chr) + (tick.dist*(0:(num.ticks-1))) - 1
      tick.labels <- sapply(tick.pos, toLabel)
      
      xplot <- ccf(chr=rep(chr.name, length(tick.pos)), x=tick.pos)$x
      y0plot <- mids(chr.name)-karyoplot$plot.params$ideogramheight/2
      if(is.null(tick.col)) {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, ...)
      } else {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, col=tick.col, ...)
      }
      graphics::text(x=xplot, y=y0plot-tick.len, labels=tick.labels, pos=1, cex=cex, offset=0.1, ...)
    
    #Minor ticks
    if(minor.ticks) {
      minor.num.ticks <- width(chr)/minor.tick.dist 
      minor.tick.pos <- start(chr) + (minor.tick.dist*(0:(minor.num.ticks-1))) - 1

      xplot <- ccf(chr=rep(chr.name, length(minor.tick.pos)), x=minor.tick.pos)$x
      y0plot <- mids(chr.name) - karyoplot$plot.params$ideogramheight/2
      if(is.null(minor.tick.col)) {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, ...)       
      } else {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, col=minor.tick.col, ...)       
      }
      
    }
  }
  
  invisible(karyoplot)
}
