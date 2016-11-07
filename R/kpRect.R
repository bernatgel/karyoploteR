#' kpRect
#' 
#' @description 
#' 
#' Plots rectangles at the specified genomic positions. 
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome context 
#' of basic plot functions
#' from R base graphics. Given a set of positions on the genome (chromosome, x0 and x1) and values 
#' (y0 and y1) for each of them, it plots rectangles going from (x0, y0) to (x1, y1). Data can be 
#' provided via a \code{GRanges} object (\code{data}), independent parameters for chr, 
#' x0, x1, y0 and y1, or a combination of both.
#' A number of parameters can be used to define exactly where and how the rectangles are drawn.
#' In addition, via the ellipsis operator (\code{...}), \code{kpRect} accepts any parameter 
#' valid for \code{rect} (e.g. \code{border}, \code{col}, ...)
#'
#' @usage kpRect(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=NULL, y1=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, ...) 
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the data. If \code{data} is present, \code{chr} will be set to \code{seqnames(data)}, \code{x0} to \code{start(data)} and x1 to \code{end(data)}. If no parameter \code{y0} is specified and \code{data} has a column named \code{y0}, this column will be used. The same for \code{y1}. (defaults to NULL)
#' @param chr    (a charecter vector) A vector of chromosome names specifying the chromosomes of the data points. If \code{data} is not NULL, \code{chr} is ignored. (defaults to NULL)
#' @param x0    (a numeric vector) A numeric vector of x left positions (in base pairs). If \code{data} is not NULL, \code{x0}. (defaults to NULL)
#' @param x1    (a numeric vector) A numeric vector of x right positions (in base pairs). If \code{data} is not NULL, \code{x1}. (defaults to NULL)
#' @param y0    (a numeric vector) A numeric vector of y bottom positions. If \code{y} is not NULL, it is used instead of any data column in \code{data}. (defaults to NULL)
#' @param y1    (a numeric vector) A numeric vector of y top positions. If \code{y} is not NULL, it is used instead of any data column in \code{data}. (defaults to NULL)
#' @param ymin    (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL)
#' @param ymax    (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpPoints}}, \code{\link{kpPlotRegions}}
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
#'   kpRect(kp, data=randomizeRegions(data.points, mask=NA), y0=0, y1=1,  r0=0, r1=0.2, border=NA, col="lightblue", data.panel=2)
#'   kpRect(kp, data=randomizeRegions(data.points, mask=NA), y0=0, y1=1,  r0=0.3, r1=0.5, border=NA, col="lightgreen", data.panel=2)
#'   kpRect(kp, data=randomizeRegions(data.points, mask=NA), y0=0, y1=1,  r0=0.6, r1=0.8, border=NA, col="purple", data.panel=2)
#'   
#' 
#'  
#' @export kpRect
#' 


kpRect <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=NULL, y1=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  pp <- prepareParameters4("kpRect", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  x0plot <- ccf(chr=pp$chr, x=pp$x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=pp$chr, x=pp$x1, data.panel=data.panel)$x
  y0plot <- ccf(chr=pp$chr, y=pp$y0, data.panel=data.panel)$y
  y1plot <- ccf(chr=pp$chr, y=pp$y1, data.panel=data.panel)$y
  
  rect(xleft=x0plot, xright=x1plot, ytop=y1plot, ybottom=y0plot, ...)
  
}
