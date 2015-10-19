


kpText <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, labels=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
 
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- kp$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- kp$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
    
  ccf <- karyoplot$coord.change.function
  
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x <- start(data) + (end(data) - start(data))/2 #Get the midpoints of the regions
    #TODO: Should we warn if regions and x are not null at the same time?
    if(is.null(y)) {
      if("value" %in% names(mcols(data))) {
        y <- data$value
      } else {
        stop("No y value specified. It is needed to provide y or a column named 'value' in data")
      }
    } 
    if(is.null(labels)) {
      if("labels" %in% names(mcols(data))) {
        labels <- data$labels
      } else {
        stop("No labels value specified. It is needed to provide labels or a column named 'labels' in data")
      }
    }
  } 
    
  #Recicle any values as needed
  chr <- recycle.first(chr, x, y, labels)
  x <- recycle.first(x, chr, y, labels)
  y <- recycle.first(y, chr, x, labels)
  labels <- recycle.first(labels, chr, x, y)
    
  
  #Scale it with ymin and ymax
  y <- (y - ymin)/(ymax - ymin)
  
  #scale y to fit in the [r0, r1] range
  y <- (y*(r1-r0))+r0
  
  xplot <- ccf(chr=chr, x=x, data.panel=data.panel)$x
  yplot <- ccf(chr=chr, y=y, data.panel=data.panel)$y
  text(x=xplot, y=yplot, labels=labels, ...)      
}
