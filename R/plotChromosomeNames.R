#internal

plotChromosomeNames <- function(karyoplot, ...) {
  chr.labels <- karyoplot$chromosomes

  x <- karyoplot$plot.params$leftmargin / 2
  y <- karyoplot$ideogram.mid(chr=chr.labels)

    #karyoplot$plot.params$xleftmargin * 0.2 + strwidth(chr.labels)/2 #Align them to the left margin, 20% into the plot

  text(x=x, y=y, labels=chr.labels, ...)
#    rect(xleft=x-strwidth(chr.labels)/2, xright=0.5+x +strwidth(chr.labels)/2, 
#         ybottom=y-strheight(chr.labels)/2, ytop=y+strheight(chr.labels)/2)
}
