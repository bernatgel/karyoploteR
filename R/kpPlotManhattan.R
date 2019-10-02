#' @name kpPlotManhattan
#' @title kpPlotManhattan
#' 
#' @description
#' Creates a manhattan plot, the ones usually seen in GWAS studies. 
#' 
#' @details 
#' Creates a manhattan plot, the ones usually seen in GWAS studies. By default,
#' it will compute the -log10 of the pvalues given and plot them as points. In
#' addition, it can plot to horizontal lines, one for the "suggestive" threshold
#' and another one for the "genomewide" significance threshold.
#' In addition, it can highlight some of the data points in a different color. 
#' Highlighted data points can be specified per name, per position in the data 
#' structure or by their position on the genome (see examples).
#' 
#' More information in the \url{https://bernatgel.github.io/karyoploter_tutorial/}{karyoploteR tutorial}.
#'  
#' @usage kpPlotManhattan(karyoplot, data=NULL, pval=NULL, points.col="2grays", 
#' points.pch=16, points.cex=1, suggestiveline = -log10(1e-05), 
#' suggestive.col="black", suggestive.lwd=1, suggestive.lty=2, 
#' genomewideline = -log10(5e-08), genomewide.col="black", genomewide.lwd=1, 
#' genomewide.lty=1, logp=TRUE, highlight=NULL, highlight.col="greenyellow",
#' ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object (or any other format accepted by \code{\link[regioneR]{toGRanges}}) with the data points to plot. If the \code{pval} parameter is NULL and  \code{data} has a column named "p" or "pval", it will be used as the pvalues to plot.
#' @param pval (numeric) The pvalues to plot. It must have the same length as data. If NULL, \code{data$p} or \code{data$pval} will be used, if present. (defaults to NULL)
#' @param points.col (colors or character) The colors used to plot the points. It can be either a vector of colors of the same length of data or a valid color specification for \code{\link{colByChr}}. (defaults to "2grays")
#' @param points.pch (numeric between 1 and 25) The symbol used to plot every point. (Defaults to 16, a filled circle)
#' @param points.cex (numeric) The size of the point symbols. (defaults to 1)
#' @param suggestiveline  (numeric) The suggestive significance threshold. The suggestive line will be plotted at this vertical position. If NULL, no line will be plotted. (defaults to -log10(1e-05))
#' @param suggestive.col (color) The color of the suggestive line. (defaults to "black")
#' @param suggestive.lwd (numeric) The width of the suggestive line (defaults to 1)
#' @param suggestive.lty (numeric) The line type of the suggestive line (defaults to 2, dashed line)
#' @param genomewideline  (numeric) The genomewide significance threshold. The genomewide line will be plotted at this vertical position. If NULL, no line will be plotted (defaults to  -log10(5e-08))
#' @param genomewide.col (color) The color of the genomewide line. (defaults to "black")
#' @param genomewide.lwd (numeric) The width of the genomewide line. (defaults to 1)
#' @param genomewide.lty (numeric) The line type of the genomewide line. (defaults to 1, solid line)
#' @param logp (logical) If TRUE, pval will be transformed using -log10(pval). (defaults to TRUE)
#' @param highlight (GRanges, character vector, logical vector or numeric vector) The points to highlight in a different color. If a GRanges (or anythng accepted by \code{\link[regioneR]{toGRanges}}) the points overlapping these regions will be highlighted. Otherwise the points will be selected with \code{data[highlight]}. If NULL no point will be highlighted. (defaults to NULL)
#' @param highlight.col The color of the highlighted points (defaults to "greenyellow")
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param ymin    (numeric) The minimum value to be plotted on the data panel. If NULL, it is set to 0. (deafults to NULL)
#' @param ymax    (numeric) The maximum value to be plotted on the data.panel. If NULL, it is set to 1. (defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions.
#'   
#' @return
#' 
#' Returns the original karyoplot object with the data computed (ymin, ymax,
#' suggestiveline, genomewideline and data (with two additional columns pval and
#' color)) stored at \code{karyoplot$latest.plot}
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPoints}}, \code{\link{colByChr}}, \code{\link{colByRegion}}
#' 
#' @examples
#' 
#' set.seed(1000)
#' 
#' #First simulate a GWAS result with a single significant peak
#' data <- createRandomRegions(nregions=20000, length.mean=1, length.sd=0, genome=filterChromosomes(getGenome("hg19")))
#' names(data) <- paste0("rs", 1:20000)
#' data$pval <- rnorm(n = 20000, mean = 0.5, sd = 0.5)
#' data$pval[data$pval<0] <- -1*data$pval[data$pval<0]
#' snps.in.peak <- which(overlapsAny(data, toGRanges("chr3:70e6-80e6")))
#' data$pval[snps.in.peak] <- runif(n = length(snps.in.peak), min=0.1, max=8)
#' data$pval <- 10^(-1*data$pval)
#'  
#' kp <- plotKaryotype("hg19", plot.type=4)
#' kp <- kpPlotManhattan(kp, data, ymax=8)
#' kpAxis(kp, ymax=8)
#' 
#' #Highlighting 
#' kp <- plotKaryotype("hg19", plot.type=4)
#' kp <- kpPlotManhattan(kp, data, ymax=8, highlight="chr3:70e6-80e6", r0=autotrack(1,4))
#' kp <- kpPlotManhattan(kp, data, ymax=8, highlight=snps.in.peak, highlight.col="orchid", r0=autotrack(2,4))
#' kp <- kpPlotManhattan(kp, data, ymax=8, highlight=names(data)[snps.in.peak], highlight.col="orange", r0=autotrack(3,4))
#' kp <- kpPlotManhattan(kp, data, ymax=8, highlight=overlapsAny(data, toGRanges("chr3:70e6-80e6")), r0=autotrack(4,4))
#' 
#' #Look and feel
#' kp <- plotKaryotype("hg19", plot.type=4)
#' kp <- kpPlotManhattan(kp, data, ymax=8, points.col="2blues", highlight="chr3:70e6-80e6", points.pch=2, points.cex=0.6, suggestive.col="red", suggestive.lwd=3, suggestive.lty=4)
#' 
#' 
#' @export kpPlotManhattan


kpPlotManhattan <- function(karyoplot, data=NULL, pval=NULL,
                            points.col="2grays", points.pch=16, points.cex=1,
                            suggestiveline = -log10(1e-05), suggestive.col="black", suggestive.lwd=1, suggestive.lty=2,
                            genomewideline = -log10(5e-08), genomewide.col="black", genomewide.lwd=1, genomewide.lty=1,
                            logp=TRUE,
                            highlight=NULL, highlight.col="greenyellow",
                            ymin=NULL, ymax=NULL, data.panel=1, 
                            r0=NULL, r1=NULL, clipping=TRUE, ...) {

  
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotManhattan: 'karyoplot' must be a valid 'KaryoPlot' object"))
  if(!methods::is(data, "GRanges")) data <- tryCatch(toGRanges(data), 
                                                     error=function(e) {stop(paste0("In kpPlotManhattan: 'data' must be a valid 'GRanges' object or a valid argument for regioneR::toGRanges"))}
                                            )

  rs <- preprocess_r0_r1(karyoplot, r0=r0, r1=r1, data.panel=data.panel)

  r0 <- rs$r0
  r1 <- rs$r1
  
  #Preprocess data
  if(is.null(pval)) {
    if("p" %in% tolower(names(mcols(data)))) {
      ncol <- which(tolower(names(mcols(data)))=="p")
      pval <- mcols(data)[,ncol]
    } else {
      if("pval" %in% tolower(names(mcols(data)))) {
        ncol <- which(tolower(names(mcols(data)))=="pval")
        pval <- mcols(data)[,ncol]
      } else {
        stop("pval was not specified.")
      }
    }
  } else {
    if(is.character(pval) && length(pval)==1) {
      pval <- data[,pval]
    }
  }
   
  if(!is.numeric(pval) || length(pval)!=length(data)) {
    stop("pval must be a numeric of the same length as data")
  }
  
  #Store pval in data (we'll return this structure at the end of the function)
  data$pval <- pval

  #Transform the pvalues
  if(logp) {
    data$y <- -1*log10(data$pval)
  } else {
    data$y <- data$pval
  }
  
  #ymax and ymin
  if(is.null(ymax)) {
    ymax <- max(max(data$y), 
                ifelse(!is.null(suggestiveline) && is.numeric(suggestiveline), suggestiveline, NULL),
                ifelse(!is.null(genomewideline) && is.numeric(genomewideline), genomewideline, NULL)
                )
  }
  if(is.null(ymin)) {
    ymin <- min(0, min(data$y), 
                ifelse(!is.null(suggestiveline) && is.numeric(suggestiveline), suggestiveline, NULL),
                ifelse(!is.null(genomewideline) && is.numeric(genomewideline), genomewideline, NULL)
            )
  }
  
  #Base Point Colors
  #If colors were specified for every point, use them 
  if(length(points.col)==length(data)) {
    data$color <- points.col
  } else {
    #Else, assume it's a color specification for colByChr
    data$color <- colByChr(data, colors = points.col)
  }
  
  #Highlights
  if(!is.null(highlight)) {
    highlight <- tryCatch(suppressWarnings(toGRanges(highlight)), error=function(e) {return(highlight)})
    if(is(highlight, "GRanges")) {
      data[overlapsAny(data, highlight)]$color <- highlight.col
    } else {
      data[highlight]$color <- highlight.col
    }
  }
  
  
  #Plot
  #TODO: remove cex, lwd, lty and col from ... if present
  
  #Plot significance lines
  if(!is.null(suggestiveline) && is.numeric(suggestiveline)) {
    kpAbline(karyoplot, h=suggestiveline, col=suggestive.col, lwd=suggestive.lwd, lty=suggestive.lty, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=TRUE, ...)
  }
  if(!is.null(genomewideline) && is.numeric(genomewideline)) {
    kpAbline(karyoplot, h=genomewideline, col=genomewide.col, lwd=genomewide.lwd, lty=genomewide.lty, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=TRUE, ...)
  }
  
  #Plot the points
  kpPoints(karyoplot, data=data, y=data$y, col=data$color, pch=points.pch, cex=points.cex, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping, ...)

  karyoplot$latest.plot <- list(funct="kpPlotManhattan", 
                                computed.values=list(ymin=ymin, ymax=ymax,
                                                     suggestiveline=suggestiveline,
                                                     genomewideline=genomewideline,
                                                     data=data))

  invisible(karyoplot)
}

