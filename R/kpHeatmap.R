#' kpHeatmap
#' 
#' @description 
#' 
#' Plots the given data as a heatmap along the genome
#' 
#' @details 
#'  
#' Given regions of the genome with a start, end and a value, draws a heatmap-like representation, with the 
#' color of the region determined by its value. It is important to note that \code{kpHeatmap} will not extend
#' the regions in any way, so if regions are not contiguous, they will appear as a series of rectangles and
#' not as a continuous plot.
#' 
#' @usage kpHeatmap(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, colors=c("blue", "white", "yellow"), ...)
#'  
#' @inheritParams kpPoints
#' @param colors    (colors) A set of color used to determine the color associated with each value. Internally, it uses \code{\link[grDevices]{colorRamp}}. (defaults to c("blue", "white", "yellow"))
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
#' 
#' set.seed(1000)
#' windows <- unlist(tileGenome(tilewidth=1000000, seqlengths = seqlengths(Hsapiens)))
#' y <- sin(x=c(1:length(windows))/10)
#' 
#' kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#' 
#' kpLines(kp, windows, y=y, r0=0.4, r1=0.6, ymin=-1, ymax=1)
#' kpAxis(kp, r0=0.4, r1=0.6, ymin=-1, ymax=1, cex=0.5)
#' 
#' kpHeatmap(kp, windows, y=y, colors = c("red", "black", "green"), r0=0, r1=0.2)
#' kpHeatmap(kp, windows, y=y, colors = c("green", "black", "red"), r0=0.2, r1=0.4)
#' 
#' #or we can provide all data into a single GRanges object
#' mcols(windows) <- data.frame(y=y)
#' 
#' kpHeatmap(kp, windows, r0=0.6, r1=0.8)
#' #non-contiguous regions appear as solitary rectangles
#' kpHeatmap(kp, sample(x = windows, 500), r0=0.8, r1=1, color=c("orange", "black", "purple", "green"))
#' 
#' 
#'@export kpHeatmap

kpHeatmap <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, colors=c("blue", "white", "yellow"), ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  ccf <- karyoplot$coord.change.function
    
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x0 <- start(data)
    x1 <- end(data)
        
    if(is.null(y)) {
      if("value" %in% names(mcols(data))) {
        y <- data$value
      } else {
        if("y" %in% names(mcols(data))) {
          y <- data$y
        } else {
          stop("No y value specified. It is needed to provide ymax or a column named 'y' in data")
        }
      }
    }
  } 
  
  if(is.null(chr)) stop("chr must be specified, either by the 'chr' parameter or by providing a 'data' object")
    
  if(is.null(ymin)) ymin <- min(y)
  if(is.null(ymax)) ymax <- max(y)
  
  if(ymin == ymax) {
    ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
    ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  }
      
  #Standardize the values to the [0,1] range to plot it with colorRamp
  y <- y - ymin
  y <- y/(ymax - ymin)
  
  
  #Create the colorRamp
  cr <- colorRamp(colors=colors)
  
 
  #Determine the plotting coordinates
  x0plot <- ccf(chr=chr, x=x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=chr, x=x1, data.panel=data.panel)$x
  yminplot <- ccf(chr=chr, y=rep_len(r0, length(chr)), data.panel=data.panel)$y
  ymaxplot <- ccf(chr=chr, y=rep_len(r1, length(chr)), data.panel=data.panel)$y
  
  rect(xleft=x0plot, xright=x1plot, ytop=ymaxplot, ybottom=yminplot, col=rgb(cr(y), max=255), border=NA, ...)
  
  invisible(karyoplot)
}
