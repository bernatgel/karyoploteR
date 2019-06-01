#' kpPlotHorizon
#' 
#' @description 
#' 
#' Plots a line joining the data points along the genome and fills the Horizon below the line.
#' 
#' @details 
#'  
#' This is a karyoploteR low-level plotting functions. Given a set of positions 
#' on the genome (chromosome and base) and a value (y) for each of them, it 
#' plots a line joining them and shades the Horizon below them. Data can be 
#' provided via a \code{GRanges} object (\code{data}), independent parameters 
#' for chr, x and y or a combination of both. A number of parameters can be used
#' to define exactly where and how the line and Horizon are drawn. In addition, 
#' via the ellipsis operator (\code{...}), \code{kpPlotHorizon} accepts any parameter 
#' valid for \code{\link[graphics]{lines}} and \code{\link[graphics]{polygon}}
#' (e.g. \code{lwd}, \code{lty}, \code{col}, \code{density}...) The lines are drawn in a per 
#' chromosome basis, so it is not possible to draw lines encompassing more than 
#' one chromosome.
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
#' set.seed(1000)
#' data.points <- sort(createRandomRegions(nregions=500, mask=NA))
#' mcols(data.points) <- data.frame(y=runif(500, min=0, max=1))
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'   kpDataBackground(kp, data.panel=1)
#'   kpDataBackground(kp, data.panel=2)
#' 
#'   kpPlotHorizon(kp, data=data.points)
#'   kpPlotHorizon(kp, data=data.points, col="lightgray", border="red", lty=2, r0=0, r1=0.5)
#'   kpPlotHorizon(kp, data=data.points, border="red", data.panel=2, r0=0, r1=0.5)
#'   kpPlotHorizon(kp, data=data.points, border="blue", data.panel=2, r0=0, r1=0.5, base.y=1)
#'
#'   kpPlotHorizon(kp, data=data.points, border="gold", data.panel=2, r0=0.5, r1=1, base.y=0.5)
#' 
#' 
#'  
#' @export kpPlotHorizon
#' 


#QUESTION: How should axis and kpHorizon relate? Should we return the values in latest plot and help creating a legend for it? )with no axis?

kpPlotHorizon <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, num.parts=3, breaks=NULL, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=NULL, r1=NULL, col="redblue6", border=NULL, clipping=TRUE, ...) {
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
      if(is.null(num.parts)) stop("In kpPlotHorizon: If breaks is NULL, num.parts cannot be NULL")
      if(!is.numeric(num.parts) || (round(num.parts)!=num.parts) || num.parts<1) stop("In kpPlotHorizon: If breaks is NULL, num.parts must be a positive integer")
      breaks <- list(pos=seq_len(num.parts-1)*ymax/num.parts,
                     neg=seq_len(num.parts-1)*ymin/num.parts)
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
  
  #kpPoints(kp, data[data$y %in% c(breaks$neg, 0, breaks$pos)], col="green", ymin=-10, ymax=10)
  
  #TODO: move to colors and export
  horizon.colors <- function(col, num.parts) {
    if(is.character(col) && length(col)==1) col <- .karyoploter.colors$horizon$schemas[[col]]
    ramp <- colorRampPalette(col)
    pal <- ramp(num.parts*2+1)
    return(list(neg=rev(pal[1:num.parts]),
                pos=pal[(num.parts+2):(2*num.parts+1)]
          ))
  }
  
  colors <- horizon.colors(col, num.parts)
  
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
