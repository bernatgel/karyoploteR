


kpBars <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, ymax=NULL, ymin=NULL, data.panel=1, bar.width=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  ccf <- karyoplot$coord.change.function
  

  
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x0 <- start(data)
    x1 <- end(data)
        
    if(is.null(ymax)) {
      if("value" %in% names(mcols(data))) {
        ymax <- data$value
      } else {
        if("ymax" %in% names(mcols(data))) {
          ymax <- data$ymax
        } else {
          stop("No y value specified. It is needed to provide y or a column named 'value' in data")
        }
      }
    }
    if(is.null(ymin)) {
      if("ymin" %in% names(mcols(data))) {
        ymin <- data$ymin
      } else {
        if(is.null(ymin)) ymin <- rep(karyoplot$plot.params[[paste0("data", data.panel, "min")]], length(ymax))
      }
    }
  } 
  
  

  x0plot <- ccf(chr=chr, x=x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=chr, x=x1, data.panel=data.panel)$x
  yminplot <- ccf(chr=chr, y=ymin, data.panel=data.panel)$y
  ymaxplot <- ccf(chr=chr, y=ymax, data.panel=data.panel)$y
  
  
  rect(xleft=x0plot, xright=x1plot, ytop=ymaxplot, ybottom=yminplot, ...)
  
}
