#' kpPlotRegions
#' 
#' @description 
#' 
#' Plots rectangles along the genome representing the regions (or intervals) specified by a \code{GRanges} object
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of karyoploteR. It takes a \code{GRanges} object and
#'  plots its content. Overlapping regions can be stacked and the number of layers for overlapping regions can be set.
#'  In contrast with the low-level functions such as \code{\link{kpRect}}, it is not possible to specify the data using 
#'  independent numeric vectors and the function only takes in \code{GRanges}.
#'
#' @usage kpPlotRegions(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col="black", border=NULL, avoid.overlapping=TRUE, num.layers=NULL, layer.margin=0.05, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the regions to plot.
# #removed as requested by the package reviewer. It can be any of the formats accepted by the \code{\link[regioneR]{toGRanges}} function from the package \href{http://bioconductor.org/packages/release/bioc/html/regioneR.html}{regioneR}.
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5 percent of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param col    (color) The background color of the regions. (defaults to black)
#' @param border    (color) The color used to draw the border of the regions. If NULL, no border is drawn. (defaults to NULL)
#' @param avoid.overlapping    (boolean) Whether overlapping regions should be drawn as stacks (TRUE) on drawing one occluding the other in a single layer (FALSE). (defaults to TRUE)
#' @param num.layers    (numeric) The number of layers the plotting space should be divided into to allow for plotting overlapping regions. The lotting region will be divided into this many pieces regardless if any overlapping regions actually exist. If NULL, the maximum number of regions overlapping a single point in the genome. (defaults to NULL)
#' @param layer.margin    (numeric) The blank space left between layers of regions. (defaults to 0.05)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#'  
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpSegments}}
#' 
#' @examples
#'  
#'  
#'  set.seed(1000)
#'  
#'  #Example 1: create 20 sets of non-overlapping random regions and plot them all. Add a coverage plot on top.
#'  kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#'  
#'  all.regs <- GRanges()
#'  
#'  nreps <- 20
#'  for(i in 1:nreps) {
#'    regs <- createRandomRegions(nregions = 100, length.mean = 10000000, length.sd = 1000000,
#'                                non.overlapping = TRUE, genome = "hg19", mask=NA)
#'    all.regs <- c(all.regs, regs)
#'    kpPlotRegions(kp, regs, r0 = (i-1)*(0.8/nreps), r1 = (i)*(0.8/nreps), col="#AAAAAA")
#'  }
#'  
#'  kpPlotCoverage(kp, all.regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
#'  kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
#'  
#'  
#'  #Example 2: Do the same with a single bigger set of possibly overlapping regions
#'  
#'  kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#'  
#'  regs <- createRandomRegions(nregions = 1000, length.mean = 10000000, length.sd = 1000000,
#'                              non.overlapping = FALSE, genome = "hg19", mask=NA)
#'                              
#'  kpPlotRegions(kp, regs, r0 = 0, r1 = 0.8, col="#AAAAAA")
#'  
#'  kpPlotCoverage(kp, regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
#'  kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
#'  
#'  
#'  
#'@export kpPlotRegions


kpPlotRegions <- function(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, 
                          autotrack=NULL, col="black", 
                          border=NULL, avoid.overlapping=TRUE, num.layers=NULL,
                          layer.margin=0.05, clipping=TRUE, ...) {
  #karyoplot
    if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(missing(data)) stop("The parameter 'data' is required")
    if(!methods::is(data, "GRanges")) stop("'data' must be a GRanges object")
  
  #if there's nothing to plot, return
  if(length(data)==0) {
    invisible(karyoplot)
  }
  
  
  #TODO: update to use toGRanges internally. But if the file is a bed file and
  # a tabix index is available, load only the regions visible in the plot.region
  # using rtracklayer::import. We might not need to check if the index is available,
  # since no overhead is expected in trying to load only a part of the bed file.
  
  if(is.null(border)) border <- col
    
  chr <- as.character(seqnames(data))
  x0 <- start(data)
  x1 <- end(data)
    
    
  #If needed, determine the y values using disjointBins to avoid region overlapping
  if(avoid.overlapping==TRUE) {
    #get the binning to avoid the overlapping
    bins <- disjointBins(data)
    if(is.null(num.layers)) num.layers <- max(bins)
    layer.height <- (1-((num.layers-1)*layer.margin))/num.layers
    y0 <- (layer.height+layer.margin) * (bins-1)
    y1 <- layer.height * (bins) + layer.margin * (bins-1)
  } else {
    #All regions in the same layer going from 0 to 1
    y0 <- 0
    y1 <- 1
  }
    
    
    
  kpRect(karyoplot=karyoplot, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=0, ymax=1, 
         r0=r0, r1=r1, autotrack=autotrack, data.panel=data.panel, col=col,
         border=border, clipping=clipping, ... )
  
  invisible(karyoplot)
}
