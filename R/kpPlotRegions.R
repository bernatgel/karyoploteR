


kpPlotRegions <- function(kp, regions, y=50, height=20, color="red", ...) {
  ccf <- kp$coord.change.function
  
  xleft <- ccf(x=start(regions))$x
  xright <- ccf(x=end(regions))$x
  ytop <- ccf(chr=as.character(seqnames(regions)), y=rep(y+height/2, length(regions)))$y
  ybottom <- ccf(chr=as.character(seqnames(regions)), y=rep(y-height/2, length(regions)))$y
  rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col=color)
  
}