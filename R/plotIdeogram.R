

plotIdeogram <- function(karyoplot, color.table=NULL, ...) {
  ccf <- karyoplot$coord.change.function

  color.table <- getIdeogramColors(color.table)
  
  if(!is.null(karyoplot$cytobands)) {
    if(length(karyoplot$cytobands)>0) { #If there are cytobands to plot, plot them

      ybottom <- ccf(as.character(seqnames(cytobands)))$y + pp$ybelowmargin
      ytop <- ccf(as.character(seqnames(cytobands)))$y + pp$ybelowmargin + pp$ideogramheight
      
      xleft <- ccf(x=start(cytobands))$x
      xright <- ccf(x=end(cytobands))$x
      
      col <- do.call(c, color.table[cytobands$gieStain])
      
      rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
      
    }
  }
    
}
