


kpPlotRegions <- function(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col=NULL, border=NULL, ...) {
  #karyoplot
    if(!hasArg(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(!hasArg(data)) stop("The parameter 'data' is required")
    
  data <- toGRanges(data)
      
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- kp$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- kp$plot.params[[paste0("data", data.panel, "max")]]
    
  ccf <- karyoplot$coord.change.function
  
  if(is.null(col)) col <- "red"
  if(is.null(border)) border <- "red"
    
  chr <- as.character(seqnames(data))
  x0 <- start(data)
  x1 <- end(data)
    
  kpRect(karyoplot=karyoplot, chr=chr, x0=x0, x1=x1, y0=0, y1=1, ymin=0, ymax=1, 
         r0=r0, r1=r1, data.panel=data.panel, col=col, border=border, ... )
}
