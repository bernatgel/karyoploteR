#' kpPlotBigWig
#' 
#' @description 
#' 
#' Plots the wiggle values in a BigWig file
#' 
#' @details 
#'  
#' \code{kpPlotBigWig} plots the data contained in a binary file format called 
#' BigWig. BigWig are used to efficiently store numeric values computed for 
#' windows covering the whole genome, ususally the coverage from an NGS 
#' experiment such as ChIP-seq. Only data required for the plotted region
#' is loaded, and when more than one chromosome is visible, it will load the
#' data for one crhomosome at a time.
#' The function accepts either a \code{\link{BigWigFile}} oject or a 
#' \code{character} with the path to a valid big wig file. The character can 
#' also be a URL to a remote server. In this case data will be loaded 
#' transparently using the \code{import} function from 
#' \code{rtracklayer}.  
#' The data is plotted using \code{\link{kpArea}} and therefore it is possible
#' to plot as a single line, a line with shaded area below or as a shaded area 
#' only adjusting the \code{col} and \code{border} parameters. 
#' 
#' 
#' @usage kpPlotBigWig(karyoplot, data, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col=NULL, border=NULL, clipping=TRUE, ...) 
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{BigWigFile} or character) The path to a bigwig file (either local or a URL to a remote file) or a \code{BigWigFile} object.
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5% of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, the minimum between 0 and the minimum value in the WHOLE GENOME will be used. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL the maximum between 0 and maximum value in the BigWigFile for the WHOLE GENOME is used. (defaults to NULL)
#' @param col  (color) The fill color of the area. If NULL the color will be assigned automatically, either a lighter version of the color used for the outer line or gray if the line color is not defined. If NA no area will be drawn. (defaults to NULL)
#' @param border  (color) The color of the line enclosing the area. If NULL the color will be assigned automatically, either a darker version of the color used for the area or black if col=NA. If NA no border will be drawn. (Defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions.
#'   
#' @return
#' 
#' Returns the original karyoplot object with the data computed (ymax and ymin values used) stored at \code{karyoplot$latest.plot}
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpArea}}, \code{\link{kpPlotBAMDensity}}
#' 
#' @examples
#' 
#' #Using the test BigWig file included in the rtracklayer package
#' rtrack_test_path <- system.file("tests", package = "rtracklayer")
#' test_bw <- file.path(rtrack_test_path, "test.bw")
#' 
#' 
#' kp <- plotKaryotype(zoom=toGRanges("chr19", 1000, 3000))
#' kp <- kpPlotBigWig(kp, data=test_bw, r0=0, r1=0.3)
#' kp <- kpPlotBigWig(kp, data=test_bw, r0=0.4, r1=0.7, border="red", lwd=2)
#' kp <- kpPlotBigWig(kp, data=test_bw, r0=0.8, r1=1, ymin=0, ymax=2, border="gold", col=NA)
#' 
#' 
#' 
#' @export kpPlotBigWig
#' @importFrom rtracklayer BigWigFile seqinfo summary
#' @importFrom S4Vectors intersect
#' 


kpPlotBigWig <- function(karyoplot, data, ymin=NULL, ymax=NULL, data.panel=1, 
                         r0=NULL, r1=NULL, autotrack=NULL, 
                         col=NULL, border=NULL, clipping=TRUE, ...) {
  
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  if(!methods::is(data, "BigWigFile") & !methods::is(data, "character")) stop("'data' must be a character or a BigWigFile object")
  
  #Prepare the access to the data
  if(methods::is(data, "character")) {
    data <- rtracklayer::BigWigFile(data)
  }
  
  
  #Check seqinfo(data) to validate at least some intersection between the genome
  #and data chromosome names
  data.chrs <- seqnames(seqinfo(data))
  if(!any(data.chrs %in% karyoplot$chromosomes)) {
    #NOTE: We can deactivate this warning and make it fail silently. Not sure 
    #      what the best option is yet
    warning("None of the chromosome names in the data file (", 
            paste0(data.chrs, collapse = ","), 
            ") matches a chromosome name in the plotted genome (", 
            paste0(karyoplot$chromosomes, collapse=","), "). Nothing will be plotted")
  }
  
  #TODO: Add a parameter to define if the ymin/ymax values should take into account: the visible regions, the visible chromosomes or the whole genome.
  # The whole genome and whole chromosomes can be computed using the statistics in the file and accessed using summary, but it seems
  # region values must be computed after loading the data.
  
  #if ymax is NULL, set it to the maximum of the file on the whole genome using
  # the summary info from BigWigFile or 0 if all values are negative
  if(is.null(ymax)) {
    ymax <- max(0, max(mcols(unlist(summary(data, type="max")))[,1]))
  }
  #if ymin is NULL, set it to the minimum of the file on the whole genome or 0 
  #if all values are above 0
  if(is.null(ymin)) {
    ymin <- min(0, min(mcols(unlist(summary(data, type="min")))[,1]))
  }
  
  
  #Load and plot the data serially for each chromosome. This will reduce the 
  # memory load when plotting on the whole genome
  
  #Filter out the plot regions (i.e. chromosomes) not known to data. This reduces work and avoids warnings about unknown seqlevels
  plot.region <- GenomeInfoDb::keepSeqlevels(karyoplot$plot.region, value = S4Vectors::intersect(seqlevels(data), seqlevels(karyoplot$plot.region)), pruning.mode = "coarse")
  
  for(i in seq_along(plot.region)) {
    #Note: remove unneeded seqlevels when calling import to avoid a warning about unknown seqlevels
    wig.data <- rtracklayer::import(data, format = "bigWig", selection=plot.region[i])
    kpArea(karyoplot, data = wig.data, y=wig.data$score, ymin=ymin, ymax=ymax,
           data.panel=data.panel, r0=r0, r1=r1, autotrack=autotrack, 
           col=col, border=border, clipping=TRUE, ...)
  }
  
  karyoplot$latest.plot <- list(funct="kpPlotBigWig", computed.values=list(ymax=ymax, ymin=ymin))
  
  
  invisible(karyoplot)
}
