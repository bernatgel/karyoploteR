

plotChromosomeNames <- function(coordChangeFunction, genome, plot.params, ...) {
  
  chr.labels <- seqlevels(genome)
  
  y <- coordChangeFunction(chr=chr.labels)$y

  x <- plot.params$xleftmargin/2

  text(x=x, y=y, labels=chr.labels)
}