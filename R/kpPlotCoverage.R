
#'@export kpPlotCoverage

kpPlotCoverage <- function(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="blue", ymin=0, ymax=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  #karyoplot
    if(!hasArg(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(!hasArg(data)) stop("The parameter 'data' is required")
        
#   #if null, get the r0 and r1
#   if(is.null(r0)) r0 <- kp$plot.params[[paste0("data", data.panel, "min")]]
#   if(is.null(r1)) r1 <- kp$plot.params[[paste0("data", data.panel, "max")]]
#     
  ccf <- karyoplot$coord.change.function
  
    
  #Compute (if needed) the coverage
  if(!is(data, "SimpleRleList")) {  #If its not a coverage object, assume it's a valid RS and compute the coverage
    data <- toGRanges(data)
    data <- coverage(data)
    #print(class(data))
  } else {#If it's already a coverage object, do nothing
    #do nothing
  }
  
  #Transform to plot
  ends <- cumsum(runLength(data))
  starts <- lapply(ends, function(x) {return(c(1, (x[-length(x)]+1)))})
  ends <- unlist(ends)  #<- This does not work in the vignette code!! Hackish way to do that follows
  starts <- starts
  #TODO: Add the first base of each chromosome to allow for chromosomes not starting at 1
  chr <- names(ends)
  value <- unlist(runValue(data))
  print(class(value))
  max.value <- max(value)
  
  print(head(chr))
  print("starts")
  print(head(unlist(starts)))
  print("ends")
  print(head(unlist(ends)))
  print("values")
  print(head(unlist(value)))
  
  if(is.null(ymax)) ymax <- max.value
  #Plot it  
   kpBars(karyoplot=karyoplot, chr=chr, x0=starts, x1=ends, y0=0, y1=value, ymin=ymin, ymax=ymax, 
          r0=r0, r1=r1, data.panel=data.panel, col=col, border=col, ... )
}
