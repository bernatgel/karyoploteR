#' kpDataBackground
#' 
#' @description 
#' 
#' Draws a solid rectangle delimiting the plotting area
#' 
#' @details 
#'  
#'  This function is used to add a background color to delimit the plotting area. 
#'  It can either delimit the whole plotting area or part of it so different data plotting
#'  regions can be seen. 
#' 
#' @usage kpDataBackground(karyoplot, r0=NULL, r1=NULL, autotrack=NULL, data.panel=1, color="gray90", clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will spàn from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5% of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param color    (color) a valid color specification
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data background will be not drawn out of the drawing area (i.e. in margins, etc) even if it overflows the visible drawing area. If FALSE, the data background representation may overflow into the margins of the plot. (defaults to TRUE)
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
#' kpAxis(kp, data.panel = 2, r0=0, r1=0.45, ymin = 0, ymax = 1, cex=0.5, 
#'        tick.pos = c(0.3, 0.5, 0.7), labels = c("-1 sd", "mean", "+1 sd"))
#' kpAxis(kp, data.panel = 2, r0=0, r1=0.45, ymin = 0, ymax = 1, cex=0.5, side=2)
#' 
#' kpDataBackground(kp, data.panel=2, r0 = 0.55, r1 = 1, color = "#EEFFEE")
#' kpAxis(kp, data.panel = 2, r0=1, r1=0.55, ymin = 0, ymax = 1, side=1, cex=0.5)
#' kpAxis(kp, data.panel = 2, r0=1, r1=0.55, ymin = 0, ymax = 1, side=2, cex=0.5)
#' 
#' #With autotrack
#' kp <- plotKaryotype("hg19", chr="chr1")
#' for(i in 1:6) {
#'   kpDataBackground(kp, autotrack = list(i, 6), color = rainbow(6)[i])
#' }
#' 
#' 
#' 
#' @export kpDataBackground

kpDataBackground <- function(karyoplot, r0=NULL, r1=NULL, autotrack=NULL, data.panel=1, color="gray90", clipping=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  #Process autotrack
  if(!is.null(autotrack) && !is.na(autotrack)) {
    rr <- processAutotrack(r0=r0, r1=r1, autotrack=autotrack)
    r0 <- rr["r0"]
    r1 <- rr["r1"]
  }
  
  xleft <- ccf(chr=as.character(seqnames(karyoplot$genome)), x=start(karyoplot$genome), data.panel=data.panel)$x
  xright <- ccf(chr=as.character(seqnames(karyoplot$genome)), x=end(karyoplot$genome), data.panel=data.panel)$x
  ytop <- ccf(chr=as.character(seqnames(karyoplot$genome)), 
              y=rep(r0, length(karyoplot$genome)), data.panel=data.panel)$y
  ybottom <- ccf(chr=as.character(seqnames(karyoplot$genome)),
                 y=rep(r1, length(karyoplot$genome)), data.panel=data.panel)$y
  
  if(karyoplot$zoom==TRUE) {
    if(clipping==TRUE) {
      dpbb <- karyoplot$getDataPanelBoundingBox(data.panel)
      graphics::clip(x1 = dpbb$xleft, x2 = dpbb$xright, y1 = dpbb$ybottom, y2=dpbb$ytop)
    }
  }
  
  graphics::rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col=color, border=FALSE, ...)
 
  invisible(karyoplot)
}