#internal

plotChromosomeNames <- function(karyoplot, ...) {
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  chr.labels <- karyoplot$chromosomes

  x <- karyoplot$plot.params$leftmargin / 2
  y <- karyoplot$ideogram.mid(chr=chr.labels)

  graphics::text(x=x, y=y, labels=chr.labels, ...)
}

