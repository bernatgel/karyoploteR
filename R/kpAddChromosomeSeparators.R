#' kpAddChromosomeSeparators
#' 
#' @description 
#' 
#' Plots between the chromosomes
#' 
#' @details 
#' 
#' Depending on the plot type it will draw vertical lines (if all chromosomes
#' are in a one line (3,4,5,7)) or horizontal lines (1,2,6)
#' 
#' By default the lines will occupy the whole chromsome extent (data.panel="all") 
#' but using the \code{data.panel} parameter it can be tuned.
#' 
#' @usage kpAddChromosomeSeparators(karyoplot, col="gray", lty=3, data.panel="all", ...) {
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param col (color) The color of the separator lines (defaults to "gray")
#' @param lty (integer) The line type of the separators (defaults to 3, dashed lines)
#' @param data.panel (data panel specification) For vertical lines, the span of the separator lines. Ignored for horizontal lines. (defaults to "all")
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}
#'  
#' @examples
#'
#' kp <- plotKaryotype(plot.type=4)
#' kpAddChromosomeSeparators(kp)
#' 
#' kp <- plotKaryotype(plot.type=5, ideogram.plotter=NULL)
#' kpAddChromosomeSeparators(kp)
#' 
#' kp <- plotKaryotype(plot.type=2)
#' kpAddChromosomeSeparators(kp, col="red")
#' 
#'  
#' @export kpAddChromosomeSeparators
#' 

kpAddChromosomeSeparators <- function(karyoplot, col="gray", lty=3, data.panel="all", ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")

  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  #If plot.type %in% 1,2,6 draw horizontal lines, if in 3,4,5,7 draw vertical lines
  ccf <- karyoplot$coord.change.function
  
  if(karyoplot$plot.type %in% c(3,4,5,7)) { #All chromosomes in one line, draw vertical separators
    y0 <- ccf(y = karyoplot$plot.params[[paste0("data", data.panel, "min")]], data.panel = data.panel)$y
    y1 <- ccf(y = karyoplot$plot.params[[paste0("data", data.panel, "max")]], data.panel = data.panel)$y
    #plot a line after the end of each chromosome except the last one
    for(nreg in seq_len(length(karyoplot$plot.region)-1)) {
      x <- ccf(chr = as.character(seqnames(karyoplot$plot.region[nreg])), x=end(karyoplot$plot.region[nreg]), data.panel=data.panel)$x
      x <- x + karyoplot$plot.params$ideogramlateralmargin/2
     
      graphics::segments(x0=x, x1=x, y0=y0, y1=y1, col = col, lty=lty, ...)
      
    }
  } else {
    #Plot horizontal lines?
    #or do nothing. I don't like them... but for consistency might be good to have it plot "something"
    
    x0 <- karyoplot$plot.params$leftmargin
    x1 <- 1 - karyoplot$plot.params$rightmargin #Could be set to the end of the shortest chromosome for each adjacent pair
    #plot a line before the each chromosome except the first one
    for(nreg in 1+seq_len(length(karyoplot$plot.region)-1)) {
      y <- ccf(chr = as.character(seqnames(karyoplot$plot.region[nreg])), y=karyoplot$plot.params$dataallmax, data.panel = "all")$y
      y <- y + karyoplot$plot.params$data1outmargin/2
      
      graphics::segments(x0=x0, x1=x1, y0=y, y1=y, col = col, lty=lty, ...)
      
    }
  }
  

  invisible(karyoplot)
}

