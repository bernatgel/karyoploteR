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
#' @usage kpAddBaseNumbers(karyoplot, tick.dist=20000000, tick.len=5, add.units=FALSE, digits=2, minor.ticks=TRUE, minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot  (karyoplot object) A valid karyoplot object created by a call to \code{\link{plotKaryotype}}
#' @param tick.dist  (numeric) The distance between the major numbered tick marks in bases (defaults to 20 milions, one major tick every 20Mb)
#' @param tick.len  (numeric) The length of the major tick marks in plot coordinates (defaults to 5)
#' @param add.units (boolean) Add the units (Mb, Kb...) to the tick labels. (Defaults to FALSE)
#' @param digits   (integer) The maximum number of digits after the decimal point in labels. (defaults to 2)
#' @param minor.ticks (boolean) Whether to add unlabeled minor ticks between the major ticks (defaults to TRUE)
#' @param minor.tick.dist   (numeric) The distance between the minor ticks in bases (defaults to 5 milions, a minor tick mark every 5Mb)
#' @param minor.tick.len   (numeric) The length of the minor tick marks in plot coordinates (defaults to 2)
#' @param cex  (numeric) The cex parameter for the major ticks label (defaults to 0.5)
#' @param tick.col   (color) If specified, the color to plot the major ticks. Otherwise the default color or, if given, the col parameter will be used. (Defaults to NULL)
#' @param minor.tick.col  (color) If specified, the color to plot the minor ticks. Otherwise the default color or, if given, the col parameter will be used. (Defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
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


kpAddBaseNumbers <- function(karyoplot, tick.dist=20000000, tick.len=5, add.units=FALSE,
                             digits=2, minor.ticks=TRUE, 
                            minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, 
                            tick.col=NULL, minor.tick.col=NULL, clipping=TRUE,  ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  
  toLabel <- function(n, add.units, digits) {
    if(add.units==TRUE) {
      units <- c("b", "Kb", "Mb")
    } else {
      units <- c("", "", "")
    }
    if(abs(n) < 1000) return(paste0(as.character(n), units[1]))
    if(abs(n) < 1000000) return(paste0(as.character(round(n/1000, digits=digits)), units[2])) #Kb
    return(paste0(as.character(round(n/1000000, digits=digits)), units[3])) #Mb
  }
  
  old.scipen <- options("scipen")
  options(scipen=999)
  on.exit(options(scipen=old.scipen), add=TRUE)

 
  #For every chromsome
  for(chr.name in karyoplot$chromosomes) {

    chr <- karyoplot$genome[chr.name]
    #Major ticks
    num.ticks <- width(chr)/tick.dist + 1
  
    tick.pos <- start(chr) + (tick.dist*(0:(num.ticks-1))) - 1 
    tick.pos[1] <- start(chr)
    
    #if zoomed in, keep only the ticks in the plot region
    if(karyoplot$zoom==TRUE) {
      if(clipping==TRUE) {
        tick.pos <- tick.pos[tick.pos >= start(karyoplot$plot.region) & tick.pos<= end(karyoplot$plot.region)]
      }
    }
    
    if(length(tick.pos)>0) {#We have to test here and cannot test on num.ticks to take the zooming into account
      tick.labels <- sapply(tick.pos, toLabel, add.units=add.units, digits=digits)
      
      
      xplot <- ccf(chr=rep(chr.name, length(tick.pos)), x=tick.pos, data.panel="ideogram")$x
      y0plot <- mids(chr.name)-karyoplot$plot.params$ideogramheight/2
      if(is.null(tick.col)) {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, ...)
      } else {
        graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len, col=tick.col, ...)
      }
      graphics::text(x=xplot, y=y0plot-tick.len, labels=tick.labels, pos=1, cex=cex, offset=0.1, ...)
    }
    #Minor ticks
    if(minor.ticks) {
      if(width(chr)>minor.tick.dist) {
        minor.num.ticks <- floor(width(chr)/minor.tick.dist)
        minor.tick.pos <- start(chr) + (minor.tick.dist*(seq_len(minor.num.ticks))) - 1
  
        #if zoomed in, keep only the ticks in the plot region
        if(karyoplot$zoom==TRUE) {
          if(clipping==TRUE) {
          minor.tick.pos <- minor.tick.pos[minor.tick.pos >= start(karyoplot$plot.region) & minor.tick.pos<= end(karyoplot$plot.region)]
          }
        }
        if(length(minor.tick.pos)>0) { 
          xplot <- ccf(chr=rep(chr.name, length(minor.tick.pos)), x=minor.tick.pos , data.panel="ideogram")$x
          y0plot <- mids(chr.name) - karyoplot$plot.params$ideogramheight/2
          if(is.null(minor.tick.col)) {
            graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, ...)       
          } else {
            graphics::segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len, col=minor.tick.col, ...)       
          }
        }
      }
    }
  }
  
  invisible(karyoplot)
}
