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
#' To do that, it will cut the \code{y} space into N parts and assign a 
#' different color to each one, flip the negative values into the positive space
#' (with negative value colors) and plot all parts in the same vertical space.
#' 
#' A more detailed explanation of horizon plots can be found at  
#' https://flowingdata.com/2015/07/02/changing-price-of-food-items-and-horizon-graphs/ 
#' and a more detailed explanation of the horizon plots implemented in 
#' karyoploteR together with an explicative animation can be found at
#' https://bernatgel.github.io/karyoploter_tutorial/ 
#' 
#' 
#' @usage kpPlotHorizon(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, 
#'                      num.parts=3, breaks=NULL, ymin=NULL, ymax=NULL,
#'                      data.panel=1, r0=0, r1=1, col="redblue6", 
#'                      border=NA, clipping=TRUE, ...)
#' 
#' @inheritParams kpPoints 
#' @param num.parts (numeric) The number of parts into which the positive and the negative y spaces will be cut. Only used if \code{breaks} is NULL. (defaults to 3)
#' @param breaks (numeric vector or list) A numeric vector of a list with two numeric vectors named "pos" and "neg". The exact break points where the y space will be cut. If NULL, the breaks will be automatically computed using num.parts. (defaults to NULL)
#' @param col  (character, color vector or list) The palette name to be used, a color vector with the color to be used to define the color gradients or a list with two color vectors named "pos" and "neg". Available palettes are: redblue6, bluepurple10 and bluegold3 (defaults to "redblue6")
#' @param border  (color) The color of the line delimiting the filled areas. If NULL the color will be assigned automatically to a darker version of the color used for the area. If NA no border will be drawn. (Defaults to NA)
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
#' kpPlotHorizon(kp, data=data.points, num.parts=1, r0=autotrack(1, 6)$r0, r1=autotrack(1, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, num.parts=2, r0=autotrack(2, 6)$r0, r1=autotrack(2, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, num.parts=3, r0=autotrack(3, 6)$r0, r1=autotrack(3, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, num.parts=5, r0=autotrack(4, 6)$r0, r1=autotrack(4, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, num.parts=9, r0=autotrack(5, 6)$r0, r1=autotrack(5, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, num.parts=15, r0=autotrack(6, 6)$r0, r1=autotrack(6, 6)$r1)
#' kpLines(kp, data=data.points, ymin=-3, ymax=3)
#' kpAbline(kp, h=0, ymin=-3, ymax=3)
#'  
#' kp <- plotKaryotype(chromosomes=c("chr1", "chr2"))
#' kpPlotHorizon(kp, data=data.points, num.parts=4, r0=autotrack(1, 6)$r0, r1=autotrack(1, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, col="redblue6", num.parts=4, r0=autotrack(2, 6)$r0, r1=autotrack(2, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, col="bluepurple10", num.parts=4, r0=autotrack(3, 6)$r0, r1=autotrack(3, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, col="bluegold3", num.parts=4, r0=autotrack(4, 6)$r0, r1=autotrack(4, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, col=c("red", "black", "green"), num.parts=4, r0=autotrack(5, 6)$r0, r1=autotrack(5, 6)$r1)
#' kpPlotHorizon(kp, data=data.points, col=c("red", "yellow"), num.parts=4, r0=autotrack(6, 6)$r0, r1=autotrack(6, 6)$r1)
#' 
#'
#' @export kpPlotHorizon
#' 


#QUESTION: How should axis and kpHorizon relate? Should we return the values in latest plot and help creating a legend for it? with no axis?

kpPlotHorizon <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, num.parts=3, breaks=NULL, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=0, r1=1, col="redblue6", border=NA, clipping=TRUE, ...) {
  #Check parameters
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotHorizon: 'karyoplot' must be a valid 'KaryoPlot' object"))
  
  #Preprocess the r's (needed to manipulate the later)
  rs <- preprocess_r0_r1(karyoplot = karyoplot, r0=r0, r1=r1, data.panel=data.panel)
  
  
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
        #If it's a numeric vector, arrange it into a valid list separating positive and negative values
        if(length(breaks)==0 || (length(breaks)>0 && all(is.numeric(breaks)))) {
          breaks <- list(pos=sort(unique(breaks[breaks>0])), neg=sort(unique(breaks[breaks<0])))
        } else {
          stop("In kpPlotHorizon: breaks must be either NULL, a numeric vector or a list with two numeric vectors called 'pos' and 'neg'")
        }
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
    colors <- horizonColors(col, max(length(breaks$pos), length(breaks$neg))+1) #num.parts)
  }
  
  #Iterate through the pos/neg regions
  for(posneg in c("pos", "neg")) {
    cols <- colors[[posneg]]
    for(nbreak in seq_len(length(breaks[[posneg]])+1)) {
      thr.1 <- c(0, breaks[[posneg]])[nbreak]
      thr.2 <- c(breaks[[posneg]], ifelse(posneg=="pos", ymax, ymin))[nbreak]
      
      d <- data
      if(posneg=="pos") {
        d$y[d$y>thr.2] <- thr.2
        d$y[d$y<thr.1] <- thr.1
        kpArea(karyoplot, d, ymin = thr.1, ymax=thr.2, col = cols[nbreak], border=border, base.y = thr.1, r0=rs$r0, r1=rs$r1)
      } else {
        d$y[d$y>thr.1] <- thr.1
        d$y[d$y<thr.2] <- thr.2
        kpArea(karyoplot, d, ymin = thr.2, ymax=thr.1, col = cols[nbreak], border=border, base.y = thr.1, r0=rs$r1, r1=rs$r0)
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
