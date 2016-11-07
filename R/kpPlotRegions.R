
#Note: y0 and y1 have no effect in kpPlotRegions

#'@export kpPlotRegions


kpPlotRegions <- function(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="black", border=NULL, avoid.overlapping=TRUE, num.layers=NULL, layer.margin=0.05, ...) {
  #karyoplot
    if(!hasArg(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(!hasArg(data)) stop("The parameter 'data' is required")
    
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
    
    
  data <- toGRanges(data)
      
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- kp$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- kp$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(border)) border <- col
    
  chr <- as.character(seqnames(data))
  x0 <- start(data)
  x1 <- end(data)
    
    
  #If needed, determine the y values using disjointBins to avoid region overlapping
  if(avoid.overlapping==TRUE) {
    #get the binning to avoid the overlapping
    bins <- disjointBins(rr)
    if(is.null(num.layers)) num.layers <- max(bins)
    layer.height <- (1-((num.layers-1)*layer.margin))/num.layers
    y0 <- (layer.height+layer.margin) * (bins-1)
    y1 <- layer.height * (bins) + layer.margin * (bins-1)
  } else {
    #All regions in the same layer going from 0 to 1
    y0 <- 0
    y1 <- 1
  }
    
    
    
  kpRect(karyoplot=karyoplot, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=0, ymax=1, 
         r0=r0, r1=r1, data.panel=data.panel, col=col, border=border, ... )
}
