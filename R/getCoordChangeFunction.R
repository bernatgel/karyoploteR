


getCoordChangeFunction <- function(genome, plot.params)
{
  #TODO: check valid params
  return(function(chr=NULL, x=NULL, y=NULL) {
    if(!is.null(y)) {
      if(is.null(chr)) {
        stop("If y is not NULL, chr must be specified too")
      }
      l <- length(chr)
      if(length(y) != l) {
        stop("If y is not NULL, it have to have the same length as chr")
      }
    }
    return(genomic2plot(chr=chr, x=x, y=y, genome=genome, plot.params=plot.params))
  })
}

getChrLowestY <- function(chr, genome, pp) {
  chr.height <- pp$ybelowmargin + pp$ideogramheight + pp$yabovemargin + pp$ydataheight
  chr.names <- seqlevels(genome)
  chrs <- c(length(chr.names):1)
  names(chrs) <- chr.names
  chr.num <- chrs[chr]
  return(pp$ybottommargin + (chr.num - 1) * chr.height)
}

#Build the function mapping genomic regions into plotting coordinates
genomic2plot <- function(chr=NULL, x=NULL, y=NULL, genome, plot.params) {
  
  pp <- plot.params
  
  genome.width <- 1 - pp$xleftmargin - pp$xrightmargin
  max.chr.len <- max(end(genome)) - min(start(genome))
    
  if(!is.null(x)) {
    x.plot <- pp$xleftmargin + (x/max.chr.len)*genome.width
  } else{
    x.plot <- NULL
  }
  
  if(is.null(chr)) {
    y.plot <- NULL 
  } else {
    chrs.y <- getChrLowestY(chr, genome, pp) 
    if(is.null(y)) { #if y is null, set y to the bottom of the chromosome
      y.plot <- chrs.y
    } else { #Return the appropiate plot.y for the given original y
      ideo.height <- pp$ybelowmargin + pp$ideogramheight + pp$yabovemargin
      datayrange <- pp$dataymax - pp$dataymin 
      yscaled <- ((y - pp$dataymin) / datayrange) * pp$ydataheight
      y.plot <- chrs.y + ideo.height + yscaled
    }
  }   
  return(list(x=x.plot, y=y.plot))
}
