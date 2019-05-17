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
#' @usage kpPlotCoverage(karyoplot, data, show.0.cov=TRUE, data.panel=1, r0=NULL, r1=NULL, col="#0e87eb", ymax=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges} or a coverage object) A GRanges object from wich the coverage will be computed or a \code{SimpleRleList} result of computing the coverage.
#' @param show.0.cov (boolean) Wether to plot a thin line representing the regions with no coverage at all. (defaults to TRUE, plot the line)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param col    (color) The background color of the regions. (defaults to "#0e87eb")
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum coverage is used. (defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#' 
#' @return
#' 
#' Returns the original karyoplot object with the data computed (max.coverage, ymax) stored in latest.plot.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotRegions}}, \code{\link{kpBars}},  \code{\link{kpPlotBAMCoverage}}, \code{\link{kpPlotDensity}}
#' 
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
#'@importFrom IRanges IntegerList
#'@importFrom GenomicRanges coverage
#'

kpPlotCoverage <- function(karyoplot, data, show.0.cov=TRUE, data.panel=1, r0=NULL, r1=NULL, col="#0e87eb", ymax=NULL, clipping=TRUE, ...) {
  #Check parameters
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  #TODO: If data is not a SimpleRleList, try to convert it to a GRanges before testing with is so other Region Set formats can be used
  
  if(!methods::is(data, "GRanges") && !methods::is(data, "SimpleRleList")) stop("'data' must be a GRanges object or a SimpleRleList")
  
  #Compute (if needed) the coverage
  #If its not a coverage object, assume it's a GRanges and compute the coverage
  if(!methods::is(data, "SimpleRleList")) { 
    #data <- toGRanges(data)
    
    #the width parameter is needed so the coverage extends to the end of the chromosomes
    data <- GenomeInfoDb::keepSeqlevels(data, karyoplot$chromosomes, pruning.mode="coarse")
    data <- GenomicRanges::coverage(data, width=karyoplot$chromosome.lengths) 
  }
  
  
  #Transform to plot
  ends <- cumsum(S4Vectors::runLength(data))
  valid.chrs <- lapply(ends, length)>0 #remove the chromosomes with no data
  ends <- ends[valid.chrs] 
  coverage.lvl <- S4Vectors::runValue(data)[valid.chrs]
  
  starts <- lapply(ends, function(x) {return(c(1, (x[-length(x)]+1)))})
  
  #convert into a GRanges
  coverage.gr <- lapply(names(ends), 
    function(chr) {
      return(toGRanges(
              data.frame(chr, starts[[chr]], ends[[chr]],
                         coverage.lvl=coverage.lvl[[chr]])
              )
             )})
  coverage.gr <- suppressWarnings(do.call(c, coverage.gr))
  
  if(show.0.cov==FALSE) {
    coverage.gr <- coverage.gr[coverage.gr$coverage.lvl!=0]
  }
  
  if(is.null(ymax)) ymax <- max(max(coverage.gr$coverage.lvl))
  
  
  kpBars(karyoplot=karyoplot, data=coverage.gr,
         y0=0, y1=coverage.gr$coverage.lvl, ymin=0, ymax=ymax,
         r0=r0, r1=r1, data.panel=data.panel,
         col="#FFAAAAAA", border=NA, clipping=clipping, ...)

  
  #NOTE: To plot with kpArea we need to build a more complex structure so
  #the actual coverage (flat tops) is drawn. I think it makes more sense to 
  #just use kpBars
  
    # cov.start <- coverage.gr
    # end(cov.start) <- start(cov.start)
    # cov.end <- coverage.gr
    # start(cov.end) <- end(cov.end)
    # cov.to.plot <- sort(c(cov.start, cov.end))
    # kpArea(karyoplot=karyoplot, data=cov.to.plot, y=cov.to.plot$coverage.lvl,
    #        ymin=0, ymax=ymax,
    #        r0=r0, r1=r1, data.panel=data.panel,
    #        col=col, border=col, clipping=clipping)
   
  karyoplot$latest.plot <- list(funct="kpPlotCoverage", computed.values=list(max.coverage=max(max(coverage.gr$coverage.lvl)),
                                                                            coverage=coverage.gr,
                                                                            ymax=ymax))
  
  invisible(karyoplot)
}
