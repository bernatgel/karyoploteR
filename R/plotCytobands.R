

plotCytobands <- function(karyoplot, color.table=NULL, add.cytobands.names=FALSE, ...) {
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
    
  color.table <- getCytobandColors(color.table)
  
  
  if(!is.null(karyoplot$cytobands)) {
    if(length(karyoplot$cytobands)>0) { #If there are cytobands to plot, plot them

      ybottom <- mids(as.character(seqnames(karyoplot$cytobands))) - pp$ideogramheight/2
      ytop <- mids(as.character(seqnames(karyoplot$cytobands))) + pp$ideogramheight/2
      
      xleft <- ccf(x=start(karyoplot$cytobands))$x
      xright <- ccf(x=end(karyoplot$cytobands))$x
      
      col <- do.call(c, color.table[karyoplot$cytobands$gieStain])
      
      rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
      
      if(add.cytobands.names) {
        plotCytobandsLabels(karyoplot=karyoplot)
      }
      
    }
  } else {
    #If no cytobands are available, plot a solid rectangle representing the chromosomes
    
    ybottom <- ccf(as.character(seqnames(karyoplot$genome)))$y + pp$ybelowmargin
    ytop <- ccf(as.character(seqnames(karyoplot$genome)))$y + pp$ybelowmargin + pp$ideogramheight
    
    xleft <- ccf(x=start(karyoplot$genome))$x
    xright <- ccf(x=end(karyoplot$genome))$x
    
    col <- "gray70"
    
    rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
    
  }
    
}
