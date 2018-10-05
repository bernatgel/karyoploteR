#' kpPlotLoess
#' 
#' @description 
#' 
#' Plot a LOESS smoothed line with confidence intervals given a list of points.
#' 
#' @details 
#'  
#' Given a set of data points (specified in any way accepted by \code{\link{kpPoints}}),
#' plot a LOESS smoothed line with optional confidence intervals. LOESS is
#' computed independently per each chromosome and data points are sorted before
#' fitting. It is possible to adjust the confidence interval with 
#' \code{conf.interval} and setting it to NULL or NA will plot no CI. It is also
#' possible to control the smoothing level with \code{span}. In addition to 
#' the standard plotting parameters, it is possible to control independently the
#' color of the fitting line and CI area and CI borders. It is also possible to
#' adjust the line type and line width of the fitting line and CI border.
#' 
#' @usage kpPlotLoess(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, conf.interval=0.95, span=0.5, data.panel=1, r0=NULL, r1=NULL, ymin=NULL, ymax=NULL, col="black", lty=1, lwd=1, ci.col="#888888AA", ci.border=NA, ci.border.lty=1, ci.border.lwd=1, clipping=TRUE, ...)
#'  
#' @inheritParams kpPoints 
#' @param conf.interval The confidence interval to plot. If a number in (0,1), the confidence interval is plotted. Else, no confidence interval is plotted. (defaults to 0.95)
#' @param span A parameter to control the smoothing level. It is passed to underlying function \code{\link[stats]{loess}}. (deafults to 0.5)
#' @param col The color of the fitting line. (defaults to "black")
#' @param lty The line type (dashed, dotted...) of the fitting line (defaults to 1, solid)
#' @param lwd The line width of the fitting line (defaults to 1)
#' @param ci.col The color of the area representing the confidence interval. If NA no CI is plotted. (defaults to #888888AA, transparent gray)
#' @param ci.border The color of line marking the border of the confidence interval. If NA, it's not plotted. (defaults to NA)
#' @param ci.border.lty The line type of line marking the border of the confidence interval. (defaults to 1, solid)
#' @param ci.border.lwd The line width of line marking the border of the confidence interval. (defaults to 1)
# @param smooth A boolean indicating if the ribbon should be smoothed (with a spline approximation) before plotting (defaults to FALSE)
#'     
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPoints}}, \code{\link{kpLines}}, \code{\link{kpPlotRibbon}}
#' 
#' @examples
#' 
#' 
#' set.seed(1000)
#' 
#' dd <- data.frame(chr="chr1", x=1:48*5e6, y=rnorm(n=48, 0.5, 0.1 ))
#' 
#' kp <- plotKaryotype(chromosomes="chr1")
#' kpPoints(kp, chr=dd$chr, x=dd$x, y=dd$y)
#' kpPlotLoess(kp, chr=dd$chr, x=dd$x, y=dd$y, col="red", conf.interval = 0.99, ci.col = "#FAAAAAAA")
#' 
#' 
#' @importFrom stats qt
#' @importFrom stats predict
#' @importFrom stats loess
#' 
#' @export kpPlotLoess
#' 
#' 


kpPlotLoess <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, conf.interval=0.95, span=0.5, data.panel=1, r0=NULL, r1=NULL, ymin=NULL, ymax=NULL, col="black", lty=1, lwd=1, ci.col="#888888AA", ci.border=NA, ci.border.lty=1, ci.border.lwd=1, clipping=TRUE, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpPlotLoess: 'karyoplot' must be a valid 'KaryoPlot' object"))

  #Use prepare parameters to normalize the input but do not rescale or filter it
  pp <- prepareParameters2("kpPlotLoess", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y,
                           ymin=0, ymax=1, r0=0, r1=1, data.panel=1, filter.data = FALSE, ...)
  
  #Prepare the arrays to accumulate the data for every chromosome to return in the karyoplot
  all.chr <- character()
  all.x <- integer()
  all.fit <- numeric()
  all.upper <- numeric()
  all.lower <- numeric()
  
  for(chr.name in karyoplot$chromosomes) {
    
    in.chr <- pp$chr==chr.name
    chr.chr <- pp$chr[in.chr]
    x.chr <- pp$x[in.chr]
    y.chr <- pp$y[in.chr]
    
    #sort the data
    ord <- order(x.chr)
    x.chr <- x.chr[ord]
    y.chr <- y.chr[ord]
    
    #compute the loess fitting with CI 
    #adapted from https://stackoverflow.com/questions/22717930/how-to-get-the-confidence-intervals-for-lowess-fit-using-r
    l.fit <- predict(loess(y.chr ~ x.chr, span=span, ...), se=TRUE)

    #Compute (and plot) the confidence interval 
    if(is.numeric(conf.interval)) {
      if(conf.interval>=1 | conf.interval<=0) {stop("Invelid value for confidence interval. Should be (0,1). Is: ", conf.interval)}
      qtile <- 1 - ((1 - conf.interval)/2)
      upper.ci <- l.fit$fit + qt(qtile, l.fit$df)*l.fit$se
      lower.ci <- l.fit$fit - qt(qtile, l.fit$df)*l.fit$se
      
      kpPlotRibbon(karyoplot=karyoplot, chr=chr.chr, x0=x.chr, x1=x.chr, y0=lower.ci, y1=upper.ci, 
                   col=ci.col, border=ci.border, lty=ci.border.lty, lwd=ci.border.lwd,
                   data.panel=data.panel, ymin=ymin, ymax=ymax, r0=r0, r1=r1,
                   clipping=clipping, ...)
    } else {
      upper.ci <- NULL
      lower.ci <- NULL
    }
    
    #And plot the fitting
    kpLines(karyoplot=karyoplot, chr=chr.chr, x=x.chr, y=l.fit$fit,
            col=col, lty=lty, lwd=lwd,
            data.panel=data.panel, ymin=ymin, ymax=ymax, r0=r0, r1=r1,
            clipping=clipping, ...)

    all.chr <- c(all.chr, chr.chr)
    all.x <- c(all.x, x.chr)
    all.fit <- c(all.fit, l.fit)
    all.upper <- c(all.upper, upper.ci)
    all.lower <- c(all.lower, upper.ci)
    
  }
    
  karyoplot$latest.plot <- list(funct="kpPlotLoess", 
                                computed.values=list(chr=all.chr,
                                                     x=all.x,
                                                     loess=all.fit, 
                                                     upper.ci=all.upper,
                                                     lower.ci=all.lower)
  )
    
  invisible(karyoplot)
}

