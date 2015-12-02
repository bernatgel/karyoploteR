


#' @export kpAxis



kpAxis <- function(karyoplot, ymin=NULL, ymax=NULL, r0=NULL, r1=NULL, side=1, numticks=3, labels=NULL, tick.pos=NULL, tick.len=NULL, label.margin=NULL, data.panel=1, ...) {
  
  
  ccf <- kp$coord.change.function
  
  if(is.null(ymin)) ymin <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(ymax)) ymax <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  if(side==1) {
    x <- start(kp$genome)
  } else {
    x <- end(kp$genome)
  }
  
  if(is.null(tick.len)) tick.len <- 0.01 * max(width(kp$genome))
  if(is.null(tick.pos)) {
    tick.pos <- (((ymax-ymin)/(numticks-1))*(0:(numticks-1)))+ymin
  } else {
    numticks <- length(tick.pos)
  }
  
  if(is.null(labels)) labels <- as.character(round(tick.pos, digits = 2))
  if(is.null(label.margin)) label.margin <-  0
  
  kpSegments(kp, chr=as.character(seqnames(kp$genome)), x0=x, x1=x, y0=ymin, y1=ymax, ymin=ymin, ymax=ymax, r0 = r0, r1=r1, data.panel=data.panel, ...)
  
 
  if(side==1) {
    kpSegments(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x0=rep(x-tick.len, each=numticks), x1=rep(x, each=numticks), y0=rep(tick.pos, length(kp$genome)), y1=rep(tick.pos, length(kp$genome)), ymin=ymin, ymax=ymax, r0 = r0, r1=r1, data.panel=data.panel, ...)
    kpText(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x=rep(x-tick.len-label.margin, each=numticks), y=rep(tick.pos, length(kp$genome)), labels = labels, ymin=ymin, ymax=ymax,  r0 = r0, r1=r1, pos=2, data.panel=data.panel, ...)  #pos=2 -> left to the given coordinate
  } else {
    kpSegments(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x0=rep(x, each=numticks), x1=rep(x+tick.len, each=numticks), y0=rep(tick.pos, length(kp$genome)), y1=rep(tick.pos, length(kp$genome)),  ymin=ymin, ymax=ymax, r0 = r0, r1=r1, ...)
    kpText(kp, chr=rep(as.character(seqnames(kp$genome)), each=numticks), x=rep(x+tick.len+label.margin, each=numticks), y=rep(tick.pos, length(kp$genome)), labels = labels,  ymin=ymin, ymax=ymax, r0 = r0, r1=r1, pos=4,  data.panel=data.panel, ...)  #pos=4 -> right to the given coordinate
  }
  
}


