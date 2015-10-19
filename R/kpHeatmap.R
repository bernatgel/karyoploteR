


kpHeatmap <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y=NULL, vmax=NULL, vmin=NULL, data.panel=1, colors=c("blue", "white", "yellow"), ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
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
  
  if(is.null(vmin)) vmin <- min(y)
  if(is.null(vmax)) vmax <- max(y)
    
  #Standardize the values to the [0,1] range to plot it with colorRamp
  y <- y - vmin
  y <- y/(vmax - vmin)
  
  
  #Create the colorRamp
  cr <- colorRamp(colors=colors)
  
  print(head(y))
  print(cr(head(y)))
  

  #Determine the plotting coordinates
  x0plot <- ccf(chr=chr, x=x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=chr, x=x1, data.panel=data.panel)$x
  yminplot <- ccf(chr=chr, y=rep_len(0, length(chr)), data.panel=data.panel)$y
  ymaxplot <- ccf(chr=chr, y=rep_len(1, length(chr)), data.panel=data.panel)$y
  
  rect(xleft=x0plot, xright=x1plot, ytop=ymaxplot, ybottom=yminplot, col=rgb(cr(y), max=255), border=NA, ...)
  
}
