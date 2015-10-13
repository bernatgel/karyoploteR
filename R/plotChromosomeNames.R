

plotChromosomeNames <- function(karyoplot, ...) {
  
  
  chr.labels <- karyoplot$chromosomes
  
  y <- karyoplot$coord.change.function(chr=chr.labels)$y

  x <- karyoplot$plot.params$xleftmargin / 2

  text(x=x, y=y, labels=chr.labels)
  
}

