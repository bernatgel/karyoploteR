#' kpBars
#' 
#' @description 
#' 
#' Plot bars along the genome
#' 
#' @details 
#'  
#' \code{kpBars} plots bars (rectangles) along the genome. It is very similar to \code{\link{kpRect}} except that if 
#' \code{y0} is missing, it's automatically set to \code{ymin} so all bars start from the base of the plotting region.
#' 
#' @usage kpBars(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y1=NULL, y0=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...)
#'  
#' @inheritParams kpRect 
#'     
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#' 
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpLines}}
#' 
#' @examples
#' 
#' 
#' set.seed(1000)
#' 
#' data <- toGRanges(data.frame(chr="chr1", start=10000000*(0:23), end=10000000*(1:24)))
#' y1 <- ((sin(start(data)) + rnorm(n=24, mean=0, sd=0.1))/5)+0.5
#' y0 <- y1 - rnorm(n=24, mean = 0, sd = 0.15)
#'  
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#' 
#' #We can specify all data values separately. If missing y0, it defaults to ymin
#' kpBars(kp, chr=as.character(seqnames(data)), x0=start(data), x1=end(data), y1=y1, col="#FFBBBB", border="#EEAAAA")
#' kpLines(kp, data=data, y=y1, col="red")
#' 
#' #or we can provide all data into a single GRanges object
#' mcols(data) <- data.frame(y0=y0, y1=y1)
#' kpBars(kp, data[data$y0>data$y1], col="orange", border="orange", data.panel=2)
#' kpBars(kp, data[data$y0<=data$y1], col="purple", border="purple", data.panel=2)
#' 
#' kpLines(kp, data, y=data$y1, data.panel=2, col="red")
#' kpLines(kp, data, y=data$y0, data.panel=2, col="blue")
#' 
#' kpAxis(kp, data.panel = 1, cex=0.8, numticks = 5, col="#777777")
#' kpAxis(kp, data.panel = 2, cex=0.8, numticks = 5, col="#777777")
#' 
#' @export kpBars


kpBars <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y1=NULL, y0=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpBars: 'karyoplot' must be a valid 'KaryoPlot' object"))
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  #If y0 is not specified with any of the valid methods, set it to the min of the data.panel
  if(is.null(y0)) {
    if(is.null(data)) {
      y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
    } else {
      if(!("y0" %in% names(mcols(data)))) {
        y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
      }
    }
  }
  
  invisible(kpRect(karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...))
  
}
