#' kpLines
#' 
#' @description 
#' 
#' Plots a line joining the data points along the genome.
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome 
#' context of basic plot functions
#' from R base graphics. Given a set of positions on the genome (chromosome and base) and a 
#' value (y) for each of them, it plots a line joining them. Data can be provided 
#' via a \code{GRanges} object (\code{data}), independent parameters for chr, x and y or a 
#' combination of both. A number of parameters can be used to define exactly where 
#' and how the lines are drawn. In addition, via the ellipsis operator (\code{...}),
#' \code{kpLines} accepts any parameter valid for \code{\link[graphics]{lines}} 
#' (e.g. \code{lwd}, \code{lty}, \code{col}, ...) The lines are drawn in a per chromosome 
#' basis, so it is not possible to draw lines encompassing more than one chromosome.
#'
#' @usage kpLines(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
#' 
#' @inheritParams kpPoints 
#' 
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
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
#'   kpPoints(kp, chr=as.character(seqnames(data.points)), 
#'            x=(start(data.points)+end(data.points))/2, y=data.points$y, pch=".",
#'            col="black", cex=1)
#' 
#'   #plotting in the data.panel=2 and using r0 and r1, ymin and ymax
#'   kpLines(kp, data=data.points, col="red", r0=0, r1=0.3, data.panel=2)
#'   kpPoints(kp, data=data.points, r0=0, r1=0.3, data.panel=2, pch=".", cex=3)
#' 
#'   kpLines(kp, data=data.points, col="blue", r0=0.4, r1=0.7, data.panel=2)
#'   kpLines(kp, data=data.points, col="blue", y=-1*(data.points$y), 
#'           ymin=-1, ymax=0, r0=0.7, r1=1, data.panel=2)
#'   #It is also possible to "flip" the data by giving an r0 > r1
#'   kpPoints(kp, data=data.points, col="red", y=(data.points$y), 
#'            r0=1, r1=0.7, data.panel=2, pch=".", cex=2)  
#' 
#' 
#'  
#' @export kpLines
#' 



kpLines <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  pp <- prepareParameters2("kpLines", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, 
                           ymin=ymin, ymax=ymax, r0=r0, r1=r1, 
                           data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  ss <- sapply(karyoplot$chromosomes, function(current.chr) {
    in.chr <- which(pp$chr==current.chr)
    xplot <- ccf(chr=pp$chr[in.chr], x=pp$x[in.chr], data.panel=data.panel)$x
    yplot <- ccf(chr=pp$chr[in.chr], y=pp$y[in.chr], data.panel=data.panel)$y
    processClipping(karyoplot=karyoplot, clipping=clipping, data.panel=data.panel)
    graphics::lines(x=xplot, y=yplot, ...)      
  })
  
  
  
  invisible(karyoplot)
}
