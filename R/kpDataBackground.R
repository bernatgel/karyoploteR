
kpDataBackground <- function(karyoplot, y0=NULL, y1=NULL, color="gray90") {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  ccf <- karyoplot$coord.change.function
  
  if(is.null(y0)) y0 <- karyoplot$plot.params$dataymin
  if(is.null(y1)) y1 <- karyoplot$plot.params$dataymax
  
  
  xleft <- ccf(x=start(karyoplot$genome))$x
  xright <- ccf(x=end(karyoplot$genome))$x
  ytop <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(y0, length(karyoplot$genome)))$y
  ybottom <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(y1, length(karyoplot$genome)))$y
  rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col=color, border=FALSE)
  
}