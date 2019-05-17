#' kpHorizon
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
#' via the ellipsis operator (\code{...}), \code{kpHorizon} accepts any parameter 
#' valid for \code{\link[graphics]{lines}} and \code{\link[graphics]{polygon}}
#' (e.g. \code{lwd}, \code{lty}, \code{col}, \code{density}...) The lines are drawn in a per 
#' chromosome basis, so it is not possible to draw lines encompassing more than 
#' one chromosome.
#' 
#   https://docs.datawatch.com/designer/tutorial/desktop/Horizon_Graph.htm
#'
#' @usage kpHorizon(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, base.y=0, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col=NULL, border=NULL, clipping=TRUE, ...)
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
#'   kpHorizon(kp, data=data.points)
#'   kpHorizon(kp, data=data.points, col="lightgray", border="red", lty=2, r0=0, r1=0.5)
#'   kpHorizon(kp, data=data.points, border="red", data.panel=2, r0=0, r1=0.5)
#'   kpHorizon(kp, data=data.points, border="blue", data.panel=2, r0=0, r1=0.5, base.y=1)
#'
#'   kpHorizon(kp, data=data.points, border="gold", data.panel=2, r0=0.5, r1=1, base.y=0.5)
#' 
#' 
#'  
#' @export kpHorizon
#' 


#QUESTION: How should axis and kpHorizon relate? Should we return the values in latest plot and help creating a legend for it? )with no axis?

kpHorizon <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, num.parts=3, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=NULL, r1=NULL, col=c("blue", "white", "red"), border=NULL, clipping=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
 
  #Normalize the parameters
  pp <- prepareParameters2("kpHorizon", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, 
                           ymin=0, ymax=1, r0=0, r1=1, data.panel=data.panel, ...)
 
  #First, find the intersections of our data with 0 and inject these positions
  #into the data to break any 0-crossing segments
  #find 0's
    shift <- function(l) {return(c(l[2:length(l)], FALSE))}
    int.0 <- which(pp$y>0 & shift(pp$y<0) | pp$y<0 & shift(pp$y>0))
      #kpPoints(kp, chr=pp$chr[int.0], x=pp$x[int.0], y=pp$y[int.0], ymin=-10, ymax=10, col="red")
    ydist <- pp$y[int.0+1] - pp$y[int.0]
    xdist <- pp$x[int.0+1] - pp$x[int.0]
    pos.0 <- pp$x[int.0] + (0-pp$y[int.0])/ydist*xdist #This is the 0's
      #kpPoints(kp, chr=pp$chr[int.0], x=pos.0, y=0, col="blue", ymin=-10, ymax=10)
  #Inject into the data
    data <- sort(c(toGRanges(pp$chr, pp$x, pp$x, y=pp$y),
                   toGRanges(pp$chr[int.0], pos.0, pos.0, y=0)
                ))
    pos.0 <- which(data$y==0)
    
    IMPORTANT: to propoerly plot this, we need to inject also the points crossing the lines
    

    #TODO: Moure a colors. Exportar? Potser si, ja posats...
    horizon.colors <- function(col, num.parts) {
      ramp <- colorRampPalette(col)
      pal <- ramp(num.parts*2+1)
      return(list(neg=rev(pal[1:num.parts]),
                  pos=pal[(num.parts+2):(2*num.parts+1)]
            ))
    }
    
    colors <- horizon.colors(col, num.parts)
    
  #Iterate through the pos/neg regions
    
    #TODO: Properly define a partition assigning a group numbering to each point and use tapply to plot
  
    last.pos <- 1
    for(curr.pos in pos.0) {
      d <- data[c(last.pos:curr.pos)]
      if(d$y[2]>0) { #if a positive range d$y[1] will be 0 usually
        range.cols <- colors$pos
      } else { #if a negative range
        range.cols <- colors$neg
        d$y <- abs(d$y)
      }
      for(i in seq_len(num.parts)) {
        y.bot <- abs((i-1)*(ymin/num.parts))
        y.top <- abs((i)*(ymin/num.parts))
        y.i <- d$y
        y.i[y.i<y.bot] <- y.bot
        y.i[y.i>y.top] <- y.top    
        kpArea(karyoplot, d, y=y.i, ymin = y.bot, ymax=y.top, col = range.cols[i], border=NA, base.y = y.bot)
      }
      last.pos <- curr.pos
    }
    
    
  #top
  
  for(i in seq_len(num.parts)) {
    y.bot <- (i-1)*(ymax/num.parts)
    y.top <- (i)*(ymax/num.parts)
    y.i <- pp$y
    y.i[y.i<y.bot] <- y.bot
    y.i[y.i>y.top] <- y.top    
    kpArea(karyoplot, chr=pp$chr, x=pp$x, y=y.i, ymin = y.bot, ymax=y.top, col = grDevices::rgb(cols(i*(1/num.parts))/255), border=NA, base.y = y.bot)
  }
  
  #bottom
  cols <- colorRamp(col[c(2,1)])
  for(i in seq_len(num.parts)) {
    y.bot <- (i-1)*(ymin/num.parts)
    y.top <- (i)*(ymin/num.parts)
    y.i <- pp$y
    y.i[y.i<y.bot] <- y.bot
    y.i[y.i>y.top] <- y.top    
    kpArea(karyoplot, chr=pp$chr, x=pp$x, y=y.i, ymin = y.bot, ymax=y.top, col = grDevices::rgb(cols(i*(1/num.parts))/255), border=NA, base.y = y.bot)
  }
  
  
     
    if(!is.na(border)) { #Plot the line only if we have a color for it
      kpLines(karyoplot, chr=pp$chr[in.chr], x= pp$x[in.chr], y=pp$y[in.chr], col=border,
              ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=clipping, ...)
    }
    
  }
  
 
  invisible(karyoplot)
}
