


kpRect <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, ymax=NULL, ymin=NULL, data.panel=1, bar.width=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  ccf <- karyoplot$coord.change.function
    
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x0 <- start(data)
    x1 <- end(data)
        
    if(is.null(ymax)) {
      if("ymax" %in% names(mcols(data))) {
        ymax <- data$ymax
      } else {
        stop("No ymax value specified. It is needed to provide ymax or a column named 'ymax' in data")
      }
    }
    if(is.null(ymin)) {
      if("ymin" %in% names(mcols(data))) {
        ymin <- data$ymin
      } else {
        stop("No ymin value specified. It is needed to provide ymin or a column named 'ymin' in data")
      }
    }
  } 
  
  #Extend ymin and ymax as necessary
  ymin <- rep_len(ymin, length(x0))
  ymax <- rep_len(ymax, length(x0))

  x0plot <- ccf(chr=chr, x=x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=chr, x=x1, data.panel=data.panel)$x
  yminplot <- ccf(chr=chr, y=ymin, data.panel=data.panel)$y
  ymaxplot <- ccf(chr=chr, y=ymax, data.panel=data.panel)$y
  
  rect(xleft=x0plot, xright=x1plot, ytop=ymaxplot, ybottom=yminplot, ...)
  
}
