#' kpAddLabels
#' 
#' @description 
#' 
#' Add labels to identify the data in the plot
#' 
#' @details 
#' 
#' Given a KaryoPlot object, plot labels on the side of the data panels to help identify the different types of data plotted
#' 
#' @usage kpAddLabels(karyoplot, labels, label.margin=0.01,  side="left", pos=NULL, offset=0, r0=NULL, r1=NULL, data.panel=1, ...) 
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param labels   (character) the text on the labels
#' @param label.margin    (numeric) the additional the margin between the labels the first base of the chromosome. In plot coordinates. Usual value might be 0.05. Can be negative. (defaults to 0.01)
#' @param side ("left" or "right") The side of the plot where to plot the labels. (defaults to "left") 
#' @param pos   (numeric) The standard graphical parameter. See \code{\link[graphics]{text}}. If NULL, pos will be selected automatically based on "side" (Defaults to NULL)
#' @param offset  (numeric) The standard graphical parameter. See \code{\link[graphics]{text}}. (Defaults to 0)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the label. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the label. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the labels are to be added. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}
#' 
#' @examples
#'
#' plot.params <- getDefaultPlotParams(plot.type=2)
#' plot.params$leftmargin <- 0.2
#' plot.params$rightmargin <- 0.2
#' 
#' #In standard whole karyotypes, labels are drawn for all chromosomes
#' 
#' kp <- plotKaryotype("hg19", chromosomes=c("chr1", "chr2"), plot.type=2, plot.params = plot.params)
#' #data panel 1
#' kpDataBackground(kp, r0=0, r1=0.5, col="#FFDDDD")
#' kpDataBackground(kp, r0=0.5, r1=1, col="#DDFFDD")
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=3, cex=0.8)
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6)
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6)
#' #data panel 2
#' kpDataBackground(kp, col="#DDDDFF", data.panel = 2)
#' kpAddLabels(kp, "BLUE", data.panel=2)
#' 
#' #Plot on the right
#' #data panel 1
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=1, cex=0.8, side="right")
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6, side="right")
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6, side="right")
#'  
#'  
#'  
#' #In karyotypes with all chromosomes in a single line, 
#' #labels are added on the first (side="left") or last (side="right") chromosome
#' 
#' kp <- plotKaryotype("hg19", chromosomes=c("chr1", "chr2", "chr3"), plot.type=3, plot.params = plot.params)
#' #data panel 1
#' kpDataBackground(kp, r0=0, r1=0.5, col="#FFDDDD")
#' kpDataBackground(kp, r0=0.5, r1=1, col="#DDFFDD")
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=3, cex=0.8)
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6)
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6)
#' #data panel 2
#' kpDataBackground(kp, col="#DDDDFF", data.panel = 2)
#' kpAddLabels(kp, "BLUE", data.panel=2)
#' 
#' #Plot on the right
#' #data panel 1
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=1, cex=0.8, side="right")
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6, side="right")
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6, side="right")
#'
#'
#'
#' #In Zoomed regions, they are placed at the correct position too
#' kp <- plotKaryotype("hg19", zoom="chr1:20000000-40000000", plot.type=2, plot.params = plot.params)
#' kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE)
#' #data panel 1
#' kpDataBackground(kp, r0=0, r1=0.5, col="#FFDDDD")
#' kpDataBackground(kp, r0=0.5, r1=1, col="#DDFFDD")
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=3, cex=0.8)
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6)
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6)
#' #data panel 2
#' kpDataBackground(kp, col="#DDDDFF", data.panel = 2)
#' kpAddLabels(kp, "BLUE", data.panel=2)
#' 
#' #Plot on the right
#' #data panel 1
#' kpAddLabels(kp, "Everything", label.margin = 0.12, srt=90, pos=1, cex=0.8, side="right")
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6, side="right")
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6, side="right")
#' 
#' 
#' 
#' 
#' @export kpAddLabels
#' 

kpAddLabels <- function(karyoplot, labels, label.margin=0.01,  side="left", pos=NULL, offset=0, r0=NULL, r1=NULL, data.panel=1, ...) {
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #Determine the position of the labels
  
  #Depending on the plot type, the labels will be printed once or multiple times
  if(karyoplot$plot.type %in% c(1,2,6)) {
    chrs <- karyoplot$chromosomes
  } 
  if(karyoplot$plot.type %in% c(3,4,5,7)) {
    if(side=="left") {
      chrs <- karyoplot$chromosomes[1]
    } else {
      chrs <- karyoplot$chromosomes[length(karyoplot$chromosomes)]
    }
  }
  
  #Compute the bounding boxes
  adj.y <- prepareParameters2("kpAddLabels", karyoplot, data=NULL, chr=chrs, 
                              x=0, y=c(0,1), ymax=1, ymin=0, r0=r0, r1=r1, 
                              data.panel=data.panel)$y
  #Compute the vertical position of the labels
  y0 <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[1], length(chrs)), data.panel = data.panel)$y
  y1 <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[2], length(chrs)), data.panel = data.panel)$y
  y <- (y0+y1)/2
  
  #Compute the horizontal position of the labels
  if(side == "left") { #if side="left", plot them on the left
    #The right end of the left-side margin minus the label.margin
    x <- rep(karyoplot$plot.params$leftmargin, length(chrs)) - label.margin
    if(is.null(pos)) pos <- 2
  } else {  #if side="right", plot them on the right
    #The last position of each chrs plus the label.margin
    x <- karyoplot$coord.change.function(chr = chrs, 
                                         x = end(karyoplot$plot.region[chrs]),
                                         y=rep(adj.y[1], length(chrs)), 
                                         data.panel = data.panel)$x + label.margin
    if(is.null(pos)) pos <- 4
  }

  graphics::text(x=x, y=y, labels=labels, pos=pos, offset=0, ...)

  invisible(karyoplot)
}

