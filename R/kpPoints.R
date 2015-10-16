


kpPoints <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, data.panel=1, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
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
  } 
    
  xplot <- ccf(chr=chr, x=x, data.panel=data.panel)$x
  yplot <- ccf(chr=chr, y=y, data.panel=data.panel)$y
  points(x=xplot, y=yplot, ...)      
}
