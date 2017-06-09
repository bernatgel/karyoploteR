#' kpCoverage
#' 
#' @description 
#' 
#' Given a \code{GRanges} object, plot the coverage along the genome. 
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of karyoploteR.
#'  It takes a \code{GRanges} object and plots it's coverage, that is, the number of regions
#'  overlapping each genomic position. The input can also be a \code{SimpleRleList} resulting
#'  from computing the coverage with \code{coverage(data)}. In contrast with the low-level 
#'  functions such as \code{\link{kpRect}}, it is not possible to specify the data using 
#'  independent numeric vectors and the function only takes in the expected object types.
#'
#' @usage kpPlotCoverage(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="#0e87eb", ymax=NULL, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object from wich the coverage will be computed or a \code{SimpleRleList} result of computing the coverage.
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param col    (color) The background color of the regions. (defaults to "#0e87eb")
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum coverage is used. (defaults to NULL)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotRegions}}, \code{\link{kpBars}}
#' @seealso \code{\link[IRanges]{coverage}}
#' 
#' @examples
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
#'  kpPlotRegions(kp, regs, r0 = 0, r1 = 0.8, col="#AAAAAA")
#'  
#'  kpPlotCoverage(kp, regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
#'  kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
#'  
#'  
#'@export kpPlotCoverage


kpPlotCoverage <- function(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="#0e87eb", ymax=NULL, ...) {
  #Check parameters
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  if(!methods::is(data, "GRanges") && !methods::is(data, "SimpleRleList")) stop("'data' must be a GRanges object or a SimpleRleList")
  
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #Compute (if needed) the coverage
  #If its not a coverage object, assume it's a GRanges and compute the coverage
  if(!methods::is(data, "SimpleRleList")) { 
    #data <- toGRanges(data)
    data <- coverage(data)
  } 
  
  #Transform to plot
  ends <- cumsum(S4Vectors::runLength(data))
  valid.chrs <- lapply(ends, length)>0 #remove the chromosomes with no data
  ends <- ends[valid.chrs] 
  coverage.lvl <- S4Vectors::runValue(data)[valid.chrs]
  
  starts <- lapply(ends, function(x) {return(c(1, (x[-length(x)]+1)))})
  
  if(is.null(ymax)) ymax <- max(max(coverage.lvl))
  
  for(chr in names(ends)) {
    kpBars(karyoplot=karyoplot, chr=chr, x0=starts[[chr]], x1=ends[[chr]], 
           y0=0, y1=coverage.lvl[[chr]], ymin=0, ymax=ymax, 
           r0=r0, r1=r1, data.panel=data.panel, col=col, border=col, ... )
  }
  
  invisible(karyoplot)
}
