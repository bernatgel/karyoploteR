#' @name kpPlotRainfall
#' @title kpPlotRainfall
#' 
#' @description 
#' 
#' Creates a rainfall plot showing the distances between features in the genome.
#' Usually used to plot the distance bewteen somatic mutations to idenify kataegis.
#' 
#' @details 
#'  
#' \code{kpPlotRainfall} plots the distances between a feature and the next one
#' in a log scale along the genome. It is usually used to plot the distance 
#' between somatic mutations in order to identify regions with an accumulation 
#' of close mutations.
#' 
#' @usage kpPlotRainfall(karyoplot, data=NULL, col=NULL, ymin=NULL, ymax=7, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the features to be plotted.
#' @param col (a color vector) The colors to use to draw the points. If the length of the vector is lower than the length of data, it will be recycled. If NULL, points will be plotted in black. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. (defaults to 7, (equivalent to 10Mb between consecutive features))
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. In particular \code{col} and \code{border} can be used to set the colors used.
#'   
#' @return
#' 
#' Returns the original karyoplot object with the data computed (distances) stored at \code{karyoplot$latest.plot}
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotDensity}}, \code{\link{kpPlotCoverage}}
#' 
#' @examples
#' 
#' set.seed(1000)
#' 
#' data <- createRandomRegions(nregions=2000)
#'  
#' kp <- plotKaryotype("hg19", plot.type=4)
#' kp <- kpPlotRainfall(kp, data)
#' kpAxis(kp, ymax=7, tick.pos=c(0:7))
#'  
#' @export kpPlotRainfall

kpPlotRainfall <- function(karyoplot, data=NULL, col=NULL, ymin=NULL, ymax=7, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotRainfall: 'karyoplot' must be a valid 'KaryoPlot' object"))
  if(!methods::is(data, "GRanges")) stop(paste0("In kpPlotRainfall: 'data' must be a valid 'GRanges' object"))

  #TODO: Convert the VCF into a GRanges if needed
  
  #Assign colors if needed
  if(is.null(col)) {
    col="#000000"
  }
  if(length(col)!=length(data)) {
    col <- rep(col, len=length(data))
  }
  
  #Sort and split by chromosome
    #first attach the color to the granges so colors are not misasigned
    mcols(data) <- data.frame(col=col, stringsAsFactors = FALSE)
    
  data <- sort(data)
  data.per.chr <- split(data, seqnames(data))
  
  distances <- list()
  #and for each chromosome, compute the distances and plot them
  for(chr in karyoplot$chromosomes) {
    if(length(data.per.chr)>0) {
      ss <- start(data.per.chr[[chr]])
      feat.dist <- ss - c(1,ss[1:length(ss)-1])
      feat.dist <- log10(feat.dist)
      distances[[chr]] <- feat.dist
      kpPoints(karyoplot, chr=chr, x=ss, y=feat.dist, col=data.per.chr[[chr]]$col, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=clipping, ...)
    }
  }
  
  karyoplot$latest.plot <- list(funct="kpPlotRainfall", computed.values=list(distances=distances))

  invisible(karyoplot)
}

