


kpHeatmap <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, colors=c("blue", "white", "yellow"), ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
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
  
}
