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
#' @usage kpAddLabels(karyoplot, labels, r0=NULL, r1=NULL, autotrack=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param labels   (character) the text on the labels
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the label. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the label. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5 percent of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param label.margin    (numeric) the additional the margin between the labels the first base of the chromosome. In plot coordinates. Usual value might be 0.05. Can be negative. (defaults to 0.01)
#' @param data.panel    (numeric) The identifier of the data panel where the labels are to be added. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param pos   (numeric) The standard graphical parameter. See \code{\link[graphics]{text}}. (Defaults to 2)
#' @param offset  (numeric) The standard graphical parameter. See \code{\link[graphics]{text}}. (Defaults to 0)
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
#' plot.params$leftmargin=0.2
#' kp <- plotKaryotype("hg19", chromosomes=c("chr1", "chr2"), plot.type=2, plot.params = plot.params)
#' #data panel 1
#' kpDataBackground(kp, r0=0, r1=0.5, col="#FFDDDD")
#' kpDataBackground(kp, r0=0.5, r1=1, col="#DDFFDD")
#' kpAddLabels(kp, "Everything", label.margin = 0.1, srt=90, pos=3, cex=0.8)
#' kpAddLabels(kp, "Red", r0=0, r1=0.5, cex=0.6)
#' kpAddLabels(kp, "Green", r0=0.5, r1=1, cex=0.6)
#' #data panel 2
#' kpDataBackground(kp, col="#DDDDFF", data.panel = 2)
#' kpAddLabels(kp, "BLUE", data.panel=2)
#'  
#' @export kpAddLabels
#' 

kpAddLabels <- function(karyoplot, labels, r0=NULL, r1=NULL, autotrack=NULL, label.margin=0.01, data.panel=1, pos=2, offset=0, ...) {
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #Determin the position of the labels
  
  #Depending on the plot type, the labels will be printed once or multiple times
  if(karyoplot$plot.type %in% c(1,2)) {
    chrs <- karyoplot$chromosomes
  } 
  if(karyoplot$plot.type %in% c(3,4,5)) {
    chrs <- karyoplot$chromosomes[1]
  }
  
  #Compute the bounding boxes
  adj.y <- prepareParameters2("kpAddLabels", karyoplot, data=NULL, chr=chrs, 
                              x=0, y=c(0,1), ymax=1, ymin=0, r0=r0, r1=r1, 
                              autotrack=autotrack, data.panel=data.panel)$y
  y0 <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[1], length(chrs)), data.panel = data.panel)$y
  y1 <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[2], length(chrs)), data.panel = data.panel)$y
  x0 <- rep(0, length(chrs))
  x1 <- rep(karyoplot$plot.params$leftmargin, length(chrs))

  x1 <- x1 - label.margin
  
  x <- x1 #So the labels are right justified
  y <- (y0+y1)/2
  
  graphics::text(x=x, y=y, labels=labels, pos=pos, offset=0, ...)

  
  invisible(karyoplot)
}

