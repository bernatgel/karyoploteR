

#internal

plotCytobandsLabels <- function(karyoplot, cytobands.names.cex=0.5, ...) {
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
      
  if(!is.null(karyoplot$cytobands)) {
    if(length(karyoplot$cytobands)>0) { #If there are cytobands to plot, plot them
  
      labels <- karyoplot$cytobands$name
      
      ylabel <- mids(as.character(seqnames(karyoplot$cytobands)))
      #ytop <- mids(as.character(seqnames(karyoplot$cytobands))) + pp$ideogramheight/2
      
      bandxleft <- ccf(x=start(karyoplot$cytobands))$x
      bandxright <- ccf(x=end(karyoplot$cytobands))$x
      
      bandmids <- (bandxleft +(bandxright-bandxleft)/2)
      
      label.width <- strwidth(labels, cex=cytobands.names.cex)
      
      do.fit <- label.width < (bandxright - bandxleft)
      #do.fit <- TRUE
      
      
      text(x=bandmids[do.fit], y=ylabel[do.fit], labels=labels[do.fit], cex=cytobands.names.cex)
      
      #rect(xleft=bandmids - label.width/2, xright=bandmids+label.width/2, ybottom=ylabel-10, ytop=ylabel+10, col=NA)      
      
    }
  }
  #If no cytobands are available, do nothing
}
