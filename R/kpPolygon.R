#' kpPolygon
#' 
#' @description 
#' 
#' Plots the the given polygons along the genome
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome context 
#' of basic plot functions from R base graphics. 
#' Given a set of positions on the genome (chromosome, base and y), it plots the polygons 
#' defined by taking these position as vertices. 
#' Data can be provided via a \code{GRanges} object (\code{data}), independent
#' parameters for chr, x and y or a combination of both. A number of parameters can be used 
#' to define exactly where and how the polygon is drawn. In addition, via the ellipsis operator 
#' (\code{...}), \code{kpPolygon} accepts any parameter valid for \code{polygon} 
#' (e.g. \code{border}, \code{density}, \code{fillOddEven}, ...)
#'
#' @usage kpPolygon(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, labels=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...)
#' 
#' @inheritParams kpPoints
#' 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpPoints}}, \code{\link{kpPlotRegions}}
#' 
#' @examples
#'  
#' set.seed(1000)
#' x <- c(1,2,5,9,13,20,15,11,7,3)*10000000
#' y <- c(0,1,0.8,0.2,0.5,0.2,1,0.3,0.1,0.2)
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'   kpDataBackground(kp, data.panel=1)
#'   kpDataBackground(kp, data.panel=2)
#' 
#'   kpPolygon(kp, chr="chr1", x=x, y=y, col="red")
#'   kpPolygon(kp, chr="chr1", x=x, y=y, col="orange", r0=0.2, r1=0.8, density=30)

#'   #use kpPolygon to draw triangles at the specified positions
#'   chr2.x <- c(1,3,7,26,48,79,120, 124, 128)*1000000
#'   for(x in chr2.x) {
#'     kpPolygon(kp, chr="chr2", x=c(x-2000000, x+2000000, x), y=c(1,1,0), r0=0, r1=0.3, col="lightblue")
#'   }
#'   
#'   
#'   kpLines(kp, data=data.points, col="blue", r0=0.4, r1=0.7, data.panel=2)
#'   kpLines(kp, data=data.points, col="blue", y=-1*(data.points$y), ymin=-1, ymax=0, r0=0.7, r1=1, data.panel=2)
#'   #It is also possible to "flip" the data by giving an r0 > r1
#'   kpPoints(kp, data=data.points, col="red", y=(data.points$y), r0=1, r1=0.7, data.panel=2, pch=".", cex=2)  
#' 

#' 
#'  
#' @export kpPolygon
#' 


kpPolygon <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  pp <- prepareParameters2("kpPolygon", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
    
  xplot <- ccf(chr=pp$chr, x=pp$x, data.panel=data.panel)$x
  yplot <- ccf(chr=pp$chr, y=pp$y, data.panel=data.panel)$y
  polygon(x=xplot, y=yplot, ...)      
}
