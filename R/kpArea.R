#' kpArea
#' 
#' @description 
#' 
#' Plots a line joining the data points along the genome and fills the area below the line.
#' 
#' @details 
#'  
#' This is a karyoploteR low-level plotting functions. Given a set of positions 
#' on the genome (chromosome and base) and a value (y) for each of them, it 
#' plots a line joining them and shades the area below them. Data can be 
#' provided via a \code{GRanges} object (\code{data}), independent parameters 
#' for chr, x and y or a combination of both. A number of parameters can be used
#' to define exactly where and how the line and area are drawn. In addition, 
#' via the ellipsis operator (\code{...}), \code{kpArea} accepts any parameter 
#' valid for \code{\link[graphics]{lines}} and \code{\link[graphics]{polygon}}
#' (e.g. \code{lwd}, \code{lty}, \code{col}, \code{density}...) The lines are drawn in a per 
#' chromosome basis, so it is not possible to draw lines encompassing more than 
#' one chromosome.
#'
#' @usage kpArea(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, base.y=0, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col=NULL, border=NULL, clipping=TRUE, ...)
#' 
#' @inheritParams kpPoints 
#' @param base.y  (numeric) The y value at which the polygon will be closed. (defaults to 0)
#' @param col  (color) The fill color of the area. A single color. If NULL the color will be assigned automatically, either a lighter version of the color used for the outer line or gray if the line color is not defined. If NA no area will be drawn. (defaults to NULL)
#' @param border  (color) The color of the line enclosing the area. A single color. If NULL the color will be assigned automatically, either a darker version of the color used for the area or black if col=NA. If NA no border will be drawn. (Defaults to NULL)
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
#'   kpArea(kp, data=data.points)
#'   kpArea(kp, data=data.points, col="lightgray", border="red", lty=2, r0=0, r1=0.5)
#'   kpArea(kp, data=data.points, border="red", data.panel=2, r0=0, r1=0.5)
#'   kpArea(kp, data=data.points, border="blue", data.panel=2, r0=0, r1=0.5, base.y=1)
#'
#'   kpArea(kp, data=data.points, border="gold", data.panel=2, r0=0.5, r1=1, base.y=0.5)
#' 
#' 
#'  
#' @export kpArea
#' 


kpArea <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, base.y=0, ymin=NULL, ymax=NULL,
                    data.panel=1, r0=NULL, r1=NULL, col=NULL, border=NULL, clipping=TRUE, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
 
  #COLORS
    #Check the length of the color parameters. it does not make sense to specify more than one in this function
    if(length(col)>1) {
      col <- col[1]
      warning("kpArea: more than one color provided in 'col'. Using only the first element ('", col, "') and ignoring the rest") 
    }
    if(length(border)>1) {
      border <- border[1]
      warning("kpArea: more than one color provided in 'border'. Using only the first element ('", border, "') and ignoring the rest") 
    }
  
    #Specify the missing colors with defaults if needed
    if(is.null(col) && (is.null(border) || is.na(border))) {
      col <- "gray70"
    }
    if(is.na(col) && is.null(border)) {
      border <- "black"
    }
    
    #And derive the missing from the other ones, if needed
    if(is.null(border) && !is.null(col) && !is.na(col)) {
      border=darker(col, amount = 100)
    }
    if(is.null(col) && !is.null(border) && !is.na(border)) {
      col=lighter(border)
    }
    
  
  #Prepare the rest of the parameters
  pp <- prepareParameters2("kpArea", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, 
                           ymin=0, ymax=1, r0=0, r1=1, data.panel=data.panel, ...)
  
  for(current.chr in karyoplot$chromosomes) {
    in.chr <- which(pp$chr==current.chr)
    
    if(!is.na(col)) { #Plot the polygon only if we have a color for it
      #prepare data for the polygon (basically add on initial and one final point at the base level)
      pol.chr <- c(current.chr, pp$chr[in.chr], current.chr)
      pol.x <- c(pp$x[in.chr][1], pp$x[in.chr], pp$x[in.chr][length(pp$x[in.chr])])
      pol.y <- c(base.y, pp$y[in.chr], base.y)
      
      kpPolygon(karyoplot, chr=pol.chr, x=pol.x, y=pol.y, col=col, border=NA, ymin=ymin, ymax=ymax,
                data.panel=data.panel, r0=r0, r1=r1, clipping=clipping, ...)
    }
    
    if(!is.na(border)) { #Plot the line only if we have a color for it
      kpLines(karyoplot, chr=pp$chr[in.chr], x= pp$x[in.chr], y=pp$y[in.chr], col=border,
              ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, clipping=clipping, ...)
    }
    
  }
  
 
  invisible(karyoplot)
}
