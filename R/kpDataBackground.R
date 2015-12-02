
#' @export kpDataBackground

kpDataBackground <- function(karyoplot, y0=NULL, y1=NULL, data.panel=1, color="gray90") {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  ccf <- karyoplot$coord.change.function
  
  if(is.null(y0)) y0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(y1)) y1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  
  xleft <- ccf(x=start(karyoplot$genome), data.panel=data.panel)$x
  xright <- ccf(x=end(karyoplot$genome), data.panel=data.panel)$x
  ytop <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(y0, length(karyoplot$genome)), data.panel=data.panel)$y
  ybottom <- ccf(chr=as.character(seqnames(karyoplot$genome)), y=rep(y1, length(karyoplot$genome)), data.panel=data.panel)$y
  rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col=color, border=FALSE)
  
}