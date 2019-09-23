#' @name kpPlotManhattan
#' @title kpPlotManhattan
#' 
#' @description 
#' 
#' Creates a rainfall plot showing the distances between features in the genome.
#' Usually used to plot the distance bewteen somatic mutations to idenify kataegis.
#' 
#' @details 
#'  
#' \code{kpPlotManhattan} plots the distances between a feature and the next one
#' in a log scale along the genome. It is usually used to plot the distance 
#' between somatic mutations in order to identify regions with an accumulation 
#' of close mutations.
#' 
#' @usage kpPlotManhattan(karyoplot, data=NULL, col=NULL, ymin=NULL, ymax=7, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
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
#' kp <- kpPlotManhattan(kp, data)
#' kpAxis(kp, ymax=7, tick.pos=c(0:7))
#'  
#' @export kpPlotManhattan

# qqman parameters
# x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
#                                                          "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
#                                                          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
#                                                          annotatePval = NULL, annotateTop = TRUE, ...
                                                         

kpPlotManhattan <- function(karyoplot, data=NULL, pval=NULL,
                            points.col="2grays", points.pch=16, points.cex=1,
                            suggestiveline = -log10(1e-05), suggestive.col="black", suggestive.lwd=1, suggestive.lty=2,
                            genomewideline = -log10(5e-08), genomewide.col="black", genomewide.lwd=1, genomewide.lty=1,
                            logp=TRUE,
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
  

  #Transform the pvalues
  if(logp) {
    pval <- -1*log10(pval)
  }
  
  
  #ymax and ymin
  if(is.null(ymax)) {
    ymax <- max(max(pval), 
                ifelse(!is.null(suggestiveline) && is.numeric(suggestiveline), suggestiveline, NULL),
                ifelse(!is.null(genomewideline) && is.numeric(genomewideline), genomewideline, NULL)
                )
  }
  if(is.null(ymin)) {
    ymin <- min(0, min(pval), 
                ifelse(!is.null(suggestiveline) && is.numeric(suggestiveline), suggestiveline, NULL),
                ifelse(!is.null(genomewideline) && is.numeric(genomewideline), genomewideline, NULL)
            )
  }
  
  #Colors
  #TODO: 
  #if it's a valid schema name, use it. If it's a color, use it.
  #Otherwhise, error
  points.col <- colByChr(data, colors = points.col)
  
  
  #Plot
  #Plot significance lines
  if(!is.null(suggestiveline) && is.numeric(suggestiveline)) {
    kpAbline(karyoplot, h=suggestiveline, col=suggestive.col, lwd=suggestive.lwd, lty=suggestive.lty, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=TRUE)
  }
  if(!is.null(genomewideline) && is.numeric(genomewideline)) {
    kpAbline(karyoplot, h=genomewideline, col=genomewide.col, lwd=genomewide.lwd, lty=genomewide.lty, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=TRUE)
  }
  
  #Plot the points
  kpPoints(karyoplot, data=data, y=pval, col=points.col, pch=points.pch, cex=points.cex, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping, ...)

  
    
  
  karyoplot$latest.plot <- list(funct="kpPlotManhattan", computed.values=list(ymin=ymin, ymax=ymax))

  invisible(karyoplot)
}

