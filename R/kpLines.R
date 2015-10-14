


kpLines <- function(karyoplot, regions=NULL, chr=NULL, x=NULL, y, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  ccf <- karyoplot$coord.change.function
  
  if(!is.null(regions)) {
    chr <- as.character(seqnames(regions))
    x <- start(regions) + (end(regions) - start(regions))/2 #Get the midpoints of the regions
    #TODO: Should we warn if regions and x are not null at the same time?
  } 
  
  ss <- sapply(seqlevels(karyoplot$genome), function(current.chr) {
    in.chr <- which(chr==current.chr)
    xplot <- ccf(chr=chr[in.chr], x=x[in.chr])$x
    yplot <- ccf(chr=chr[in.chr], y=y[in.chr])$y
    lines(x=xplot, y=yplot)      
  })  
}
