


kpAbline <- function(karyoplot, chr=NULL, h=NULL, v=NULL, ymin=NULL, ymax=NULL, data.panel=1,  r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
 
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  ccf <- karyoplot$coord.change.function
  
  if(is.null(chr)) { #if chr is not specified, plot the line in all chromosomes
    chr <- as.character(seqnames(karyoplot$genome)) 
  } 
    
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  
  if(!is.null(h)) {
    x0 <- start(karyoplot$genome[seqnames(karyoplot$genome) == chr])
    x1 <- end(karyoplot$genome[seqnames(karyoplot$genome) == chr])
    
    kpSegments(karyoplot=karyoplot, chr=chr, x0=x0, x1=x1, y0=h, y1=h, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)  
  }    
   
  if(!is.null(v)) {
    y0 <- ymin 
    y1 <- ymax
    
    kpSegments(karyoplot=karyoplot, chr=chr, x0=v, x1=v, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)  
  }
        
}
