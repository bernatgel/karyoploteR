#' kpSegments
#' 
#' @description 
#' 
#' Plots segments at the specified genomic positions. 
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome context 
#' of basic plot functions from R base graphics. 
#' Given a set of positions on the genome (chromosome, x0 and x1) and values 
#' (y0 and y1) for each of them, it plots segments going from (x0, y0) to (x1, y1). Data can be 
#' provided via a \code{GRanges} object (\code{data}), independent parameters for chr, 
#' x0, x1, y0 and y1, or a combination of both.
#' A number of parameters can be used to define exactly where and how the segments are drawn.
#' In addition, via the ellipsis operator (\code{...}), \code{kpSegments} accepts any parameter 
#' valid for \code{segments} (e.g. \code{lwd}, \code{lty}, \code{col}, ...)
#'
#' @usage kpSegments(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=NULL, y1=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, ...) 
#' 
#' @inheritParams kpRect 
#' 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpPoints}}, \code{\link{kpPlotRegions}}
#' 
#' @examples
#'  
#' set.seed(1000)
#' data.points <- sort(createRandomRegions(nregions=500, length.mean=2000000, mask=NA))
#' y <- runif(500, min=0, max=0.8)
#' mcols(data.points) <- data.frame(y0=y, y1=y+0.2)
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'   kpDataBackground(kp, data.panel=1)
#'   kpDataBackground(kp, data.panel=2)
#' 
#'   kpRect(kp, data=data.points, col="black")
#'   kpSegments(kp, data=data.points, col="white")
#'   
#'   kpSegments(kp, data=data.points, y0=0, y1=1,  r0=0.2, r1=0.8, col="lightblue", data.panel=2)
#'   kpSegments(kp, data=data.points, y0=0, y1=1,  r0=0.8, r1=0.2, col="lightgreen", data.panel=2)
#'   

#' 
#'  
#' @export kpSegments
#' 



kpSegments <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL,  ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  pp <- prepareParameters4("kpSegments", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  x0plot <- ccf(chr=pp$chr, x=pp$x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=pp$chr, x=pp$x1, data.panel=data.panel)$x
  y0plot <- ccf(chr=pp$chr, y=pp$y0, data.panel=data.panel)$y
  y1plot <- ccf(chr=pp$chr, y=pp$y1, data.panel=data.panel)$y
  
  segments(x0=x0plot, x1=x1plot, y0=y0plot, y1=y1plot, ...)
}
