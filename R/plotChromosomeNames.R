

plotChromosomeNames <- function(karyoplot, ...) {
  
  chr.labels <- karyoplot$chromosomes
  
  y <- karyoplot$coord.change.function(chr=chr.labels)$y + karyoplot$plot.params$ybelowmargin + karyoplot$plot.params$ideogramheight/2 + karyoplot$plot.params$ydataheight/2

  x <- karyoplot$plot.params$xleftmargin / 2
    #karyoplot$plot.params$xleftmargin * 0.2 + strwidth(chr.labels)/2 #Align them to the left margin, 20% into the plot

  text(x=x, y=y, labels=chr.labels)
#   rect(xleft=x-strwidth(chr.labels)/2, xright=0.5+x +strwidth(chr.labels)/2, 
#        ybottom=y-strheight(chr.labels)/2, ytop=y+strheight(chr.labels)/2)
}
