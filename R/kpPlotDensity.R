#' kpPlotDensity
#' 
#' @description 
#' 
#' Plots the density of features along the genome
#' 
#' @details 
#'  
#' \code{kpPlotDensity} plots the density of a set of features represented by a 
#' \code{GRanges} object along the genome. It creates a non-overlapping tiling 
#' of the genome and computes the number of features per window. It's possible 
#' to specify the window size.
#' 
#' @usage kpPlotDensity(karyoplot, data=NULL, window.size=1e6, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object from which the density will be computed.
#' @param window.size (numeric) The size of the windows for wich the density is computed. (Defaults to 1e6, one megabase windows)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum density is used. (defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. In particular \code{col} and \code{border} can be used to set the colors used.
#'   
#' @return
#' 
#' Returns the original karyoplot object with the data computed (windows and density) stored at \code{karyoplot$latest.plot}
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotRibbon}}, \code{\link{kpPlotCoverage}}
#' 
#' @examples
#' 
#' set.seed(1000)
#' 
#' data <- createRandomRegions(nregions=20000)
#'  
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes="chr1")
#' 
#' kp <- kpPlotDensity(kp, data)
#' kpAxis(kp, ymin = 0, ymax=kp$latest.plot$computed.values$max.density)
#' 
#' kp <- kpPlotDensity(kp, data, data.panel=2, col="#CCCCFF",  ymax=20, lwd=2)
#' kpAxis(kp, ymin = 0, ymax=20, data.panel=2)
#'
#' kp <- kpLines(kp, data=kp$latest.plot$computed.values$windows, y=kp$latest.plot$computed.values$density, col="black", r0=0.5, r1=1, data.panel=2, ymax=20)
#' 
#' @export kpPlotDensity


kpPlotDensity <- function(karyoplot, data=NULL, window.size=1e6, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL,  clipping=TRUE, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotDensity: 'karyoplot' must be a valid 'KaryoPlot' object"))
  if(!methods::is(data, "GRanges")) stop(paste0("In kpPlotDensity: 'data' must be a valid 'GRanges' object"))
 
  #TODO: Add more checks - window.size numeric...
  
  #create bins all along the genome
  windows <- tileGenome(
          stats::setNames(karyoplot$chromosome.lengths, karyoplot$chromosomes),
          tilewidth = window.size, cut.last.tile.in.chrom = TRUE)
  
  dens <- countOverlaps(windows, data)
  
  if(is.null(ymax)) {
    ymax <- max(dens)
  }
  
  
  karyoplot <- kpPlotRibbon(karyoplot, data = windows, y0=0, y1=dens, ymin=ymin,
                            ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, 
                            clipping=clipping, ...)

  karyoplot$latest.plot <- list(funct="kpPlotDensity", 
                                computed.values=list(density=dens, 
                                                     windows=windows,
                                                     max.density=max(dens))
                                )

  invisible(karyoplot)
}

