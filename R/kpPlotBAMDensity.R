#' kpPlotBAMDensity
#' 
#' @description 
#' 
#' Plots the density of features along the genome
#' 
#' @details 
#'  
#' \code{kpPlotBAMDensity} plots the read density of a BAM file. It does not
#' plot the coverage but the read density as the number of reads overlapping 
#' a every window. It uses \code{\link{Rsamtools}} to efficiently access the
#' BAM file. The BAM file must be indexed.
#' 
#' @usage kpPlotBAMDensity(karyoplot, data=NULL, window.size=1e6, normalize=FALSE, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col="gray80", border=NA,  clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{BamFile} or character) The path to a bam file (must be indexed) or a \code{BamFile} object.
#' @param window.size (numeric) The size of the windows for wich the density is computed. (Defaults to 1e6, one megabase windows)
#' @param normalize (boolean) Specifies if the density values should be normalized by the total number of mapped reads in the bam file. (Defaults to FALSE)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum density is used. (defaults to NULL)
#' @param col  (color) The background color to plot. If NULL, it will be a lighter version of 'border' or 'black' if border is null. (Defaults to "gray80")
#' @param border (color) The color to use to plot the borders of the bars. If NULL, it will be a darker version of 'col'. If NA, no border will be plotted. (Defaults to NULL)
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
#' library(pasillaBamSubset) #A package with 2 example bam files
#' un1.bam.file <- untreated1_chr4() # get the name of the first bam
#' un3.bam.file <- untreated3_chr4() #and the name of the second
#' 
#' window.size <- 1e4 #compute the density with 10kb windows
#' 
#' kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
#' kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
#' kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, r0=0.5, r1=1, ymax=50000, col="darkorange")
#' kp <- kpPlotBAMDensity(kp, data = un3.bam.file, window.size = window.size, r0=0.5, r1=0, ymax=50000, col="darkorchid") #using r0>r1 we can flip the plot
#' kpAxis(kp, ymin=0, ymax=50000, r0=0.5, r1=1, labels = c("0", "25K", "50K"))
#' kpAxis(kp, ymin=0, ymax=50000, r0=0.5, r1=0, labels = c("0", "25K", "50K"))
#' 
#' kpText(kp, chr = "chr4", x=7e5, y=0.85, labels = paste0("Untreated 1 (reads per ", window.size, " bases)"))
#' kpText(kp, chr = "chr4", x=7e5, y=0.15, labels = paste0("Untreated 3 (reads per ", window.size, " bases)"))
#' 
#' 
#' 
#' #Or normalizing by the number of mapped reads
#' kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
#' kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
#' kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, normalize=TRUE, r0=0.5, r1=1, ymax=0.2, col="darkorange")
#' kp <- kpPlotBAMDensity(kp, data = un3.bam.file, window.size = window.size, normalize=TRUE, r0=0.5, r1=0, ymax=0.2, col="darkorchid") #using r0>r1 we can flip the plot
#' 
#' 
#' @export kpPlotBAMDensity
#' @importFrom Rsamtools BamFile countBam idxstatsBam ScanBamParam


kpPlotBAMDensity <- function(karyoplot, data=NULL, window.size=1e6, normalize=FALSE, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col="gray80", border=NA, clipping=TRUE,...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotBAMDensity: 'karyoplot' must be a valid 'KaryoPlot' object"))
  if(is.character("data")) {
    data <- Rsamtools::BamFile(file = data)
  }
  if(!methods::is(data, "BamFile")) stop(paste0("In kpPlotBAMDensity: 'data' must be a character or a 'BamFile' object."))
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  

  #use tileGenome to create windows only on the visible part of the genome, that
  #is in the karyoplot$plot.region
  plot.region.lengths <- setNames(width(karyoplot$plot.region), as.character(seqnames(karyoplot$plot.region)))
  windows <- tileGenome(plot.region.lengths, tilewidth = window.size, cut.last.tile.in.chrom = TRUE)
  seqinfo(windows) <- Seqinfo(seqnames=seqlevels(windows)) #remove the seqlength info from seqinfo to avoid a potential out-of-bounds warning when shifting the windows
  
  #Now, move the windows start(karyoplot$plor.region) bases to the right. 
  #It's only necessary when zoomed or with chromosomes not starting at position 1
  windows <- shift(windows, shift=start(karyoplot$plot.region[seqnames(windows)])-1)
  
  
  #Count the number of read in the bam overlapping each window
  dens <- Rsamtools::countBam(data, param=Rsamtools::ScanBamParam(which = windows))$records
  
  if(normalize==TRUE) {
    total.reads <- sum(Rsamtools::idxstatsBam(file=data)$mapped)
    dens <- dens/total.reads
  }
  
  if(is.null(ymax)) {
    ymax <- max(dens)
  }
  
  #Specify the missing colors if possible
  if(is.null(border) & !is.null(col) & !is.na(col)) {
    border=darker(col, amount = 100)
  }
  if(is.null(col) & !is.null(border) & !is.na(col)) {
    col=lighter(border)
  }
  if(is.na(col) & is.null(border)) {
    border <- "black"
  }
  
  karyoplot <- kpBars(karyoplot, data=windows, y0=0, y1=dens,ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, border=border, col=col, clipping=clipping, ...)
  
  karyoplot$latest.plot <- list(funct="kpPlotBAMDensity", computed.values=list(density=dens, windows=windows, max.density=max(dens)))

  invisible(karyoplot)
}

