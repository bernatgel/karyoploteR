#' kpAxis
#' 
#' @description 
#' 
#' Plot axis at the sides of the data panels
#' 
#' @details 
#'  
#' \code{kpAxis} plots axis at the sides of the data panels. It is possible to control the number of ticks and their labels,
#'  the placement of the plots and whether they span the whole data panel or just part of it. To do that they use the same placement 
#'  parameters used by other karyoploteR functions (\code{r0} and \code{r1}). This function does not  Axis are always plotted for all chromosomes.
#' 
#' @usage kpAxis(karyoplot, ymin=NULL, ymax=NULL, r0=NULL, r1=NULL, side=1, numticks=3, labels=NULL, tick.pos=NULL, tick.len=NULL, label.margin=NULL, data.panel=1, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param ymin    (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL)
#' @param ymax    (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param side  (numeric) In which side of the data panel should the axis be plotted. 1 - plot it on the right of the data panel. 2 - Plot it on the left. (defualts to 1)
#' @param num.ticks    (numeric) the number of ticks (and labels) of the axis. If \code{tick.pos} is present, it takes precedence over \code{num.ticks} and \code{num.ticks} is ignored. (defaults to 3)
#' @param labels    (character) the labels to be placed next to the ticks. If the number of labels is lower than the number of  tickes, the labels will be reused. If NULL, the numeric values of the ticks will be used. (defaults to NULL)
#' @param tick.pos    (numeric) the places in the axis where a tick should be drawn. If present, \code{num.ticks} is ignored. If NULL, ticks are placed equidistant. (defaults to NULL)
#' @param tick.len    (numeric) the length of the ticks to be drawn measured in base pairs. If NULL, tick length is 0.01 times the length in bases of the longest chromosome. (defaults to NULL)
#' @param label.margin    (numeric) the additional the margin between the labels and ticks. Can be negative. If NULL, the default margin is used. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#'
#'    
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpDataBackground}}, \code{\link{kpAbline}}
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
#' @export kpAxis



kpAxis <- function(karyoplot, ymin=NULL, ymax=NULL, r0=NULL, r1=NULL, side=1, numticks=3, labels=NULL, tick.pos=NULL, tick.len=NULL, label.margin=NULL, data.panel=1, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- kp$coord.change.function
  
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(side==1) {
    x <- start(kp$genome)
  } else {
    x <- end(kp$genome)
  }
  
  if(is.null(tick.len)) tick.len <- 0.01 * max(width(kp$genome))
  if(is.null(tick.pos)) {
    tick.pos <- (((ymax-ymin)/(numticks-1))*(0:(numticks-1)))+ymin
  } else {
    numticks <- length(tick.pos)
  }
  
  if(is.null(labels)) labels <- as.character(round(tick.pos, digits = 2))
  if(is.null(label.margin)) label.margin <-  0
  
  kpSegments(kp, chr=as.character(seqnames(kp$genome)), x0=x, x1=x, y0=ymin, y1=ymax, ymin=ymin, ymax=ymax, r0 = r0, r1=r1, data.panel=data.panel, ...)
  
 
  if(side==1) {
    kpSegments(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x0=rep(x-tick.len, each=numticks), x1=rep(x, each=numticks), y0=rep(tick.pos, length(kp$genome)), y1=rep(tick.pos, length(kp$genome)), ymin=ymin, ymax=ymax, r0 = r0, r1=r1, data.panel=data.panel, ...)
    kpText(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x=rep(x-tick.len-label.margin, each=numticks), y=rep(tick.pos, length(kp$genome)), labels = labels, ymin=ymin, ymax=ymax,  r0 = r0, r1=r1, pos=2, data.panel=data.panel, ...)  #pos=2 -> left to the given coordinate
  } else {
    kpSegments(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x0=rep(x, each=numticks), x1=rep(x+tick.len, each=numticks), y0=rep(tick.pos, length(kp$genome)), y1=rep(tick.pos, length(kp$genome)),  ymin=ymin, ymax=ymax, r0 = r0, r1=r1, data.panel=data.panel, ...)
    kpText(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x=rep(x+tick.len+label.margin, each=numticks), y=rep(tick.pos, length(kp$genome)), labels = labels,  ymin=ymin, ymax=ymax, r0 = r0, r1=r1, pos=4,  data.panel=data.panel, ...)  #pos=4 -> right to the given coordinate
  }
  
  invisible(karyoplot)
}


