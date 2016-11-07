

#'@export plotCytobands

plotCytobands <- function(karyoplot, color.table=NULL, add.cytobands.names=FALSE, add.base.numbers=FALSE, ...) {
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
    
  color.table <- getCytobandColors(color.table)
  
  #If there are no cytobands, create a single "fake" cytoband to represent the whole chromosome
  if(!is.null(karyoplot$cytobands) && length(karyoplot$cytobands)>0) { 
    cyto <- karyoplot$cytobands  
  } else {
    cyto <- karyoplot$genome
    mcols(cyto) <- data.frame(name=seqnames(cyto), gieStain="gpos50", stringsAsFactors=FALSE)  
  }
  
  #And plot them  
  
  ybottom <- mids(as.character(seqnames(cyto))) - pp$ideogramheight/2
  ytop <- mids(as.character(seqnames(cyto))) + pp$ideogramheight/2
    
  xleft <- ccf(x=start(cyto))$x
  xright <- ccf(x=end(cyto))$x
    
  col <- do.call(c, color.table[cyto$gieStain])
    
  rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
    
  if(add.cytobands.names) {
    plotCytobandsLabels(karyoplot=karyoplot, ...)
  }
    
  if(add.base.numbers) {
    plotBaseNumbers(karyoplot, ...)
  }
}
