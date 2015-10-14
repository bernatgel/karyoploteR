

plotCytobands <- function(karyoplot, color.table=NULL, ...) {
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  
  color.table <- getCytobandColors(color.table)
  
  
  
  if(!is.null(karyoplot$cytobands)) {
    if(length(karyoplot$cytobands)>0) { #If there are cytobands to plot, plot them

      ybottom <- ccf(as.character(seqnames(karyoplot$cytobands)))$y + pp$ybelowmargin
      ytop <- ccf(as.character(seqnames(karyoplot$cytobands)))$y + pp$ybelowmargin + pp$ideogramheight
      
      xleft <- ccf(x=start(karyoplot$cytobands))$x
      xright <- ccf(x=end(karyoplot$cytobands))$x
      
      col <- do.call(c, color.table[karyoplot$cytobands$gieStain])
      
      rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
      
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
