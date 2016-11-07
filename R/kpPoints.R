#' kpPoints
#' 
#' @description 
#' 
#' Plots data points along the genome.
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome context 
#' of basic plot functions
#' from R base graphics. Given a set of positions on the genome (chromosome and base) and a 
#' value (y) for each of them, it plots the set of points representing them. Data can be provided 
#' via a \code{GRanges} object (\code{data}), independent parameters for chr, x and y or a 
#' combination of both. A number of parameters can be used to define exactly where 
#' and how the points are drawn. In addition, via the ellipsis operator (\code{...}), \code{kpPoints}
#' accepts any parameter valid for \code{points} (e.g. \code{pch}, \code{cex}, \code{col}, ...)
#'
#' @usage kpPoints(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the data. If \code{data} is present, \code{chr} will be set to \code{seqnames(data)} and \code{x} to the midpoints of the rages in data. If no parameter \code{y} is specified and \code{data} has a column named \code{y} or \code{value} this column will be used to define the \code{y} value of each data point. (defaults to NULL)
#' @param chr    (a charecter vector) A vector of chromosome names specifying the chromosomes of the data points. If \code{data} is not NULL, \code{chr} is ignored. (defaults to NULL)
#' @param x    (a numeric vector) A numeric vector with the positions (in base pairs) of the data points in the chromosomes. If \code{data} is not NULL, \code{x} is ignored. (defaults to NULL)
#' @param y    (a numeric vector) A numeric vector with the values of the data points. If \code{y} is not NULL, it is used instead of any data column in \code{data}. (defaults to NULL)
#' @param ymin    (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL)
#' @param ymax    (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpText}}, \code{\link{kpPlotRegions}}
#' 
#' @examples
#'  
#' set.seed(1000)
#' data.points <- sort(createRandomRegions(nregions=500, mask=NA))
#' mcols(data.points) <- data.frame(y=runif(500, min=0, max=1))
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'   kpDataBackground(kp, data.panel=1)
#'   kpDataBackground(kp, data.panel=2)
#' 
#'   kpLines(kp, data=data.points, col="red")
#' 
#'   #Three ways of specifying the exact same data.points
#'   kpPoints(kp, data=data.points)
#'   kpPoints(kp, data=data.points, y=data.points$y, pch=16, col="#CCCCFF", cex=0.6)
#'   kpPoints(kp, chr=as.character(seqnames(data.points)), x=(start(data.points)+end(data.points))/2, y=data.points$y, pch=".", col="black", cex=1)
#' 
#'   #plotting in the data.panel=2 and using r0 and r1, ymin and ymax
#'   kpLines(kp, data=data.points, col="red", r0=0, r1=0.3, data.panel=2)
#'   kpPoints(kp, data=data.points, r0=0, r1=0.3, data.panel=2, pch=".", cex=3)
#' 
#'   kpLines(kp, data=data.points, col="blue", r0=0.4, r1=0.7, data.panel=2)
#'   kpLines(kp, data=data.points, col="blue", y=-1*(data.points$y), ymin=-1, ymax=0, r0=0.7, r1=1, data.panel=2)
#'   #It is also possible to "flip" the data by giving an r0 > r1
#'   kpPoints(kp, data=data.points, col="red", y=(data.points$y), r0=1, r1=0.7, data.panel=2, pch=".", cex=2)  
#' 

#' 
#'  
#' @export kpPoints
#' 


kpPoints <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  pp <- prepareParameters2("kpPoints", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  xplot <- ccf(chr=pp$chr, x=pp$x, data.panel=data.panel)$x
  yplot <- ccf(chr=pp$chr, y=pp$y, data.panel=data.panel)$y
  points(x=xplot, y=yplot, ...)   
  
  invisible(karyoplot)
}
