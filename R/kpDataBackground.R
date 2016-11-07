#' kpDataBackground
#' 
#' @description 
#' 
#' Draws a solid rectangle delimiting the plotting area
#' 
#' @details 
#'  
#'  This function is used to add a background color to delimit the plotting area. It can either delimit the
#'   whole plotting area or part of it so different data plotting regions can be seen. 
#' 
#' @usage kpDataBackground(karyoplot, r0=NULL, r1=NULL, data.panel=1, color="gray90", ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param color    (color) a valid color specification
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#'
#'    
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpAxis}}
#' 
#' @examples
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'
#' #Prepare data panel 1
#' kpDataBackground(kp, data.panel=1)
#' kpAxis(kp, data.panel = 1)
#' kpAxis(kp, data.panel = 1, ymin = 0, ymax=10, numticks = 11, side = 2, cex = 0.4, col="red")
#'
#' #Prepare data panel 2
#' #Data panel 2 is conceptually split into two parts and the second part is "inverted"
#' kpDataBackground(kp, data.panel=2, r0 = 0, r1 = 0.45, color = "#EEEEFF")
#' kpAxis(kp, data.panel = 2, r0=0, r1=0.45, ymin = 0, ymax = 1, cex=0.5, tick.pos = c(0.3, 0.5, 0.7), labels = c("-1 sd", "mean", "+1 sd"))
#' kpAxis(kp, data.panel = 2, r0=0, r1=0.45, ymin = 0, ymax = 1, cex=0.5, side=2)
#' 
#' kpDataBackground(kp, data.panel=2, r0 = 0.55, r1 = 1, color = "#EEFFEE")
#' kpAxis(kp, data.panel = 2, r0=1, r1=0.55, ymin = 0, ymax = 1, side=1, cex=0.5)
#' kpAxis(kp, data.panel = 2, r0=1, r1=0.55, ymin = 0, ymax = 1, side=2, cex=0.5)
#' 
#' 
#' @export kpDataBackground

kpDataBackground <- function(karyoplot, r0=NULL, r1=NULL, data.panel=1, color="gray90", ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  
  xleft <- ccf(x=start(karyoplot$genome), data.panel=data.panel)$x
  xright <- ccf(x=end(karyoplot$genome), data.panel=data.panel)$x
  ytop <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(r0, length(karyoplot$genome)), data.panel=data.panel)$y
  ybottom <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(r1, length(karyoplot$genome)), data.panel=data.panel)$y
  rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col=color, border=FALSE, ...)
 
  invisible(karyoplot)
}