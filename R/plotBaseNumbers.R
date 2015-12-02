

#internal

plotBaseNumbers <- function(karyoplot, tick.dist=20000000, tick.len=5, minor.ticks=TRUE, minor.tick.dist=5000000, minor.tick.len=2,  cex=0.5, ...) {
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  
  toLabel <- function(n) {
    if(abs(n) < 1000) return(as.character(n))
    if(abs(n) < 1000000) return(paste0(as.character(round(n/10)/100), "")) #Kb
    return(paste0(as.character(round(n/10000)/100), "")) #Mb
  }
  
  old.scipen <- options("scipen")
  options(scipen=999)
  on.exit(options(scipen=old.scipen), add=TRUE)
  
  #For every chromsome
  for(chr.name in as.character(seqnames(karyoplot$genome))) {
    chr <- karyoplot$genome[as.character(seqnames(karyoplot$genome)) == chr.name]
    #Major ticks
      num.ticks <- width(chr)/tick.dist 
      tick.pos <- start(chr) + (tick.dist*(0:(num.ticks-1))) - 1
      tick.labels <- sapply(tick.pos, toLabel)
      
      xplot <- ccf(chr=chr.name, x=tick.pos)$x
      y0plot <- mids(chr.name)-kp$plot.params$ideogramheight/2
      segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-tick.len)
      text(x=xplot, y=y0plot-tick.len, labels=tick.labels, pos=1, cex=cex, offset=0.1)
    
    #Minor ticks
    if(minor.ticks) {
      minor.num.ticks <- width(chr)/minor.tick.dist 
      minor.tick.pos <- start(chr) + (minor.tick.dist*(0:(minor.num.ticks-1))) - 1

      xplot <- ccf(chr=chr.name, x=minor.tick.pos)$x
      y0plot <- mids(chr.name)-kp$plot.params$ideogramheight/2
      segments(x0=xplot, x1=xplot, y0=y0plot, y1=y0plot-minor.tick.len)     
    }
  }

}
