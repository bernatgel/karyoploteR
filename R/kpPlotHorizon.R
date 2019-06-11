#' kpPlotHorizon
#' 
#' @description 
#' 
#' Plot a horizon plot, an area-like plot where different value levels are
#' plotted in different colors.
#' 
#' @details 
#'  
#' kpPlotHorizon will create a horizon plot, a plot usually used in time series,
#' that will represent a wiggle value (coverage, methylation, expression...) 
#' that we'd usually plot with a line or area in a fraction of the vertical
#' space used by these function (usually in 1/6, but that's configurable).
#' To do that, 
#' 
#   https://docs.datawatch.com/designer/tutorial/desktop/Horizon_Graph.htm
#'
#' @usage kpPlotHorizon(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, base.y=0, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col=NULL, border=NULL, clipping=TRUE, ...)
#' 
#' @inheritParams kpPoints 
#' @param base.y  (numeric) The y value at wich the polygon will be closed. (defaults to 0)
#' @param col  (color) The fill color of the Horizon. If NULL the color will be assigned automatically, either a lighter version of the color used for the outer line or gray if the line color is not defined. If NA no Horizon will be drawn. (defaults to NULL)
#' @param border  (color) The color of the line enclosing the Horizon. If NULL the color will be assigned automatically, either a darker version of the color used for the Horizon or black if col=NA. If NA no border will be drawn. (Defaults to NULL)
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpText}}, \code{\link{kpPlotRibbon}}
#' 
#' @examples
#'
#'
#' data.points <- toGRanges(data.frame(chr=c("chr1", "chr1"), start=c(1, 100), end=c(1, 100), y=c(-2, 2)))
#' kp <- plotKaryotype(zoom=toGRanges("chr1:1-100"))
#' kpLines(kp, data.points, r0=0, r1=0.5, ymin=-2, ymax=2)
#' kpAbline(kp, h=c(-2,-1,0,1,2), col="#AAAAAA", r0=0, r1=0.5, ymin=-2, ymax=2)
#' kpAxis(kp, r0=0, r1=0.5, ymin=-2, ymax=2)
#' kpPlotHorizon(kp, data.points, r0=0.55, r1=1)
#'
#'
#' data.points <- toGRanges(data.frame(chr="chr1", start=1:100, end=1:100))
#' data.points$y <- sin(start(data.points)/2)
#' kp <- plotKaryotype(zoom=toGRanges("chr1:1-100"))
#' kpLines(kp, data.points, r0=0, r1=0.5, ymin=-1, ymax=1)
#' kpAbline(kp, h=c(-1,0,1), col="#AAAAAA", r0=0, r1=0.5, ymin=-1, ymax=1)
#' kpAxis(kp, r0=0, r1=0.5, ymin=-1, ymax=1)
#' kpPlotHorizon(kp, data.points, r0=0.55, r1=1)
#'
#'
#' set.seed(1000)
#' data.points <- sort(createRandomRegions(nregions=500, mask=NA))
#' mcols(data.points) <- data.frame(y=runif(500, min=-3, max=3))
#' 
#' kp <- plotKaryotype(chromosomes=c("chr1", "chr2"))
#' kpPlotHorizon(kp, data=data.points)
#'
#' kp <- plotKaryotype(chromosomes=c("chr1", "chr2"))
#' kpPlotHorizon(kp, data=data.points, num.breaks=2, r0=0.01, r1=0.2)
#' kpPlotHorizon(kp, data=data.points, num.breaks=3, r0=0.21, r1=0.4)
#' kpPlotHorizon(kp, data=data.points, num.breaks=5, r0=0.42, r1=0.6)
#' kpPlotHorizon(kp, data=data.points, num.breaks=9, r0=0.61, r1=0.8)
#' kpPlotHorizon(kp, data=data.points, num.breaks=15, r0=0.81, r1=1)
#' kpLines(kp, data=data.points, ymin=-3, ymax=3)
#'  
#' @export kpPlotHorizon
#' 


#QUESTION: How should axis and kpHorizon relate? Should we return the values in latest plot and help creating a legend for it? )with no axis?

kpPlotHorizon <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, num.breaks=3, breaks=NULL, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=0, r1=1, col="redblue6", border=NULL, clipping=TRUE, ...) {
  #Check parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotHorizon: 'karyoplot' must be a valid 'KaryoPlot' object"))
  
  #Normalize the parameters
  pp <- prepareParameters2("kpPlotHorizon", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, 
                           ymin=0, ymax=1, r0=0, r1=1, data.panel=data.panel, ...)
  #And build a GRanges with the normalized data for easier management later on
  data <- toGRanges(pp$chr, pp$x, pp$x, y=pp$y)
  
  #Define undefined parameters
    #ymin and ymax
    if(is.null(ymin)) {
      ymin <- min(0, min(data$y))
    }
    if(is.null(ymax)) {
      ymax <- max(0, max(data$y))
    }
  
    #breaks
    if(is.null(breaks)) {
      if(is.null(num.breaks)) stop("In kpPlotHorizon: If breaks is NULL, num.parts cannot be NULL")
      if(!is.numeric(num.breaks) || (round(num.breaks)!=num.breaks) || num.breaks<1) stop("In kpPlotHorizon: If breaks is NULL, num.parts must be a positive integer")
      breaks <- list(pos=seq_len(num.breaks-1)*ymax/num.breaks,
                     neg=seq_len(num.breaks-1)*ymin/num.breaks)
    } else {
      if(!is.list(breaks) || !setequal(names(breaks), c("pos", "neg")) || !all(lapply(breaks, methods::is, "numeric"))) {
        stop("In kpPlotHorizon: breaks must be either NULL or a list with two numeric vectors called 'pos' and 'neg'")
      }
    }
  

  #Find the intersections of the data lines with the breaks and 0 and 
  #inject them into the data
  for(thr in c(breaks$neg, 0, breaks$pos)) {
    isecs <- findIntersections(data, thr)
    data <- c(data, isecs)
  }
  data <- sort(data)
  
  if(is.list(col) && c("pos", "neg" %in% names(col))) {
    if(all(is.color(col$neg)) && all(is.color(col$pos))) {
      colors <- col
    } else {
      stop("in kpPlotHorizon: col$neg and col$pos must be valid colors. ")
    }
  } else {
    colors <- horizonColors(col, num.breaks)
  }
  
  #Iterate through the pos/neg regions
  for(posneg in c("pos", "neg")) {
    cols <- colors[[posneg]]
    for(nbreak in seq_len(length(breaks$pos)+1)) {
      thr.1 <- c(0, breaks[[posneg]])[nbreak]
      thr.2 <- c(breaks[[posneg]], ifelse(posneg=="pos", ymax, ymin))[nbreak]
      
      d <- data
      if(posneg=="pos") {
        d$y[d$y>thr.2] <- thr.2
        d$y[d$y<thr.1] <- thr.1
        kpArea(karyoplot, d, ymin = thr.1, ymax=thr.2, col = cols[nbreak], border=NA, base.y = thr.1, r0=r0, r1=r1)
      } else {
        d$y[d$y>thr.1] <- thr.1
        d$y[d$y<thr.2] <- thr.2
        kpArea(karyoplot, d, ymin = thr.2, ymax=thr.1, col = cols[nbreak], border=NA, base.y = thr.1, r0=r1, r1=r0)
      }
    }
  }
  
 
  karyoplot$latest.plot <- list(funct="kpPlotHorizon", 
                                computed.values=list(breaks=breaks, 
                                                     colors=colors
                                                     )
  )
  
 
  invisible(karyoplot)
}
