

prepareParameters2 <- function(function.name, karyoplot, data, chr, x, y, ymax, ymin, r0, r1, data.panel, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop(paste0("In ", function.name, ": 'karyoplot' must be a valid 'KaryoPlot' object"))
    
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- kp$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- kp$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(!is.null(data)) {
    chr <- as.character(seqnames(data))
    x <- start(data) + (end(data) - start(data))/2 #Get the midpoints of the regions
    
    if(is.null(y)) {
      if("y" %in% names(mcols(data))) {
        y <- data$y
      } else {
        if("value" %in% names(mcols(data))) {
          y <- data$value
        } else {
          stop("No y value specified. Parameter y or a column named 'y' or 'value' in data must be provided")  
        }
      }
    } 
  } 
  
  if(is.null(chr)) stop("chr must be specified, either by the 'chr' parameter or by providing a 'data' object")
  
  #Recicle any values as needed
  chr <- recycle.first(chr, x, y)
  x <- recycle.first(x, chr, y)
  y <- recycle.first(y, chr, x)
  
  #Scale it with ymin and ymax
  y <- (y - ymin)/(ymax - ymin)
  
  
  #scale y to fit in the [r0, r1] range
  y <- (y*(r1-r0))+r0
    
  return(list(chr=chr, x=x, y=y))
  
}
