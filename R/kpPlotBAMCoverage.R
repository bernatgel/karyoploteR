#' kpPlotBAMCoverage
#' 
#' @description 
#' 
#' Plots the coverage of a BAM file along the genome
#' 
#' @details 
#'  
#' \code{kpPlotBAMCoverage} plots the read coverage of a BAM file, that is, the 
#' number of reads overlapping each position. It uses the 
#' \code{\link{bamsignals}} package to efficiently access the BAM file.
#' The BAM file must be indexed. This function is only recommended when 
#' plotting small parts of the genome. For larger plots consider using 
#' \code{\link{kpPlotBAMDensity}}.
#'  
#' 
#' @note Since the plotting the exact coverage for large
#' regions of the genome may be unfeasable, it includes a safety mechanism
#' causing it to raise a warning and do nothing if the region is larger than 
#' a threshold specified by \code{max.valid.region.size}.
#' 
#' @usage kpPlotBAMCoverage(karyoplot, data=NULL, max.valid.region.size=1e6, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col=NULL, border=NA, clipping=TRUE,...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a character) The path to a bam file (must be indexed).
#' @param max.valid.region.size    (numeric) If the length of plotted region exceeds this number, nothing will be plotted. It's a safety mechanism to from excessive memory usage. (Defaults to 1e6, 1 milion bases)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5 percent of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum density is used. (defaults to NULL)
#' @param col  (color) The fill color to plot. If NULL the color will be assigned automatically, either a lighter version of the color used for the outer line or gray if the line color is not defined. If NA no area will be drawn. (defaults to NULL)
#' @param border (color) The color to use to plot the borders of the bars. If NULL, it will be a darker version of 'col'. If NA, no border will be plotted. (Defaults to NA)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. In particular \code{col} and \code{border} can be used to set the colors used.
#'   
#' @return
#' 
#' Returns the original karyoplot object with the data computed (max.coverage) stored at \code{karyoplot$latest.plot}
#' 
#' @seealso \code{\link{kpPlotBAMDensity}}, \code{\link{kpPlotCoverage}}
#' 
#' @examples
#' 
#' library(pasillaBamSubset) #A package with 2 example bam files
#' un1.bam.file <- untreated1_chr4() # get the name of the first bam
#' un3.bam.file <- untreated3_chr4() #and the name of the second
#' 
#' kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
#' kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
#' kp <- kpPlotBAMCoverage(kp, data = un1.bam.file) #Warning and does not plot. region too large.
#' kp <- kpPlotBAMCoverage(kp, data = un1.bam.file, max.valid.region.size=2000000)
#'
#' #Use zoom to plot a smaller region to see the coverage with more detail
#' kp <- plotKaryotype(genome="dm6", zoom=toGRanges("chr4", 340000, 350000))
#' kp <- kpAddBaseNumbers(kp, tick.dist = 1e3)
#' kp <- kpPlotBAMCoverage(kp, data = un1.bam.file)
#' 
#' 
#' #Change the colors and borders and compare  two bams
#' kp <- plotKaryotype(genome="dm6", zoom=toGRanges("chr4", 340000, 350000))
#' kp <- kpAddBaseNumbers(kp, tick.dist = 1e3)
#' kp <- kpPlotBAMCoverage(kp, data = un1.bam.file, r0=0.5, r1=1, border="orange")
#' kp <- kpPlotBAMCoverage(kp, data = un3.bam.file, r0=0.5, r1=0, border="darkgreen") #r1 < r0 will flip the plot
#' kpAbline(kp, h=0.5, col="darkgray")
#' 
#' 
#' 
#' @export kpPlotBAMCoverage
#' @importFrom bamsignals bamCoverage


kpPlotBAMCoverage <- function(karyoplot, data=NULL, max.valid.region.size=1e6,
                              ymin=NULL, ymax=NULL, data.panel=1,
                              r0=NULL, r1=NULL, autotrack=NULL, 
                              col=NULL, border=NA, clipping=TRUE,...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotBAMCoverage: 'karyoplot' must be a valid 'KaryoPlot' object"))

  if(!is.character(data)) stop(paste0("In kpPlotBAMCoverage: 'data' must be a character indicating the path to an indexed BAM file."))
  if(!file.exists(data)) stop(paste0("In kpPlotBAMCoverage: File ", data, " does not exists"))
  if(!file.exists(paste0(data, ".bai"))) stop(paste0("In kpPlotBAMCoverage: Index file ", paste0(data, ".bai"), " does not exists. The BAM file must be indexed."))

  if(sum(width(karyoplot$plot.region))>max.valid.region.size) {
    warning("In kpPlotBAMCoverage: Skipping BAM coverage plot. 
The genomic region in your plot is larger than the maximum
valid size. Plotting coverage in a very large region would 
probably result in very large memory usage and potentially 
crash your R session while creating not very informative 
plots. For larger regions, kpPlotBAMDensity is recommended.
If you want to plot the coverage in a large region,
please supply a larger value for 'max.valid.region.size'
parameter.")
    invisible(karyoplot)
  }
  
  #Everything look correct. Start plotting.
  
  #Get the coverage
  bam.cov <- bamsignals::bamCoverage(data, karyoplot$plot.region, verbose=FALSE)
  #bam.cov has a numeric vector for each entry in plot.region
  
  max.cov <- max(unlist(lapply(bam.cov, max)))
  if(is.null(ymax)) {
    ymax <- max.cov
  }
  
  for(i in seq_len(length(karyoplot$plot.region))) {
    r <- karyoplot$plot.region[i]
    kpArea(karyoplot, chr=seqnames(r), x=seq.int(start(r), end(r)), 
           y=bam.cov[i], ymin=ymin, ymax=ymax, data.panel=data.panel, 
           r0=r0, r1=r1, autotrack=autotrack, border=border, col=col, clipping=clipping, ...)
  }
  
  karyoplot$latest.plot <- list(funct="kpPlotBAMCoverage", computed.values=list(max.coverage=max.cov))

  invisible(karyoplot)
}

