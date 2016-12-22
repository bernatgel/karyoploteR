#' @keywords internal

#Return the coordinate change function mapping genomic to plot coordinates


getCoordChangeFunctions <- function(plot.type, genome, plot.params)
{
 
  if(plot.type == 1) {
    genomic2plot <- genomic2plot_1HorizDataAboveIdeogram
    ideoMid <- getIdeogramMidY_1HorizDataAboveIdeogram
    chrHeight <- getChrHeight_1HorizDataAboveIdeogram
  }
  if(plot.type == 2) {
    genomic2plot <- genomic2plot_2HorizDataAboveAndBelowIdeogram
    ideoMid <- getIdeogramMidY_2HorizDataAboveAndBelowIdeogram
    chrHeight <- getChrHeight_2HorizDataAboveAndBelowIdeogram
  }
  if(plot.type == 4) {
    genomic2plot <- genomic2plot_4VerticalDataAboveAndBelowIdeogram
    ideoMid <- getIdeogramMidY_4VerticalDataAboveAndBelowIdeogram
    chrHeight <- getChrHeight_4VerticalDataAboveAndBelowIdeogram
  }
  
  #TODO: check valid params
  return(list(
    coorChangeFunction=function(chr=NULL, x=NULL, y=NULL, data.panel=NULL) {
    if(!is.null(y)) {
      if(is.null(chr)) {
        stop("If y is not NULL, chr must be specified too")
      }
      l <- length(chr)
      if(length(y) != l) {
        stop("If y is not NULL, it have to have the same length as chr")
      }
    }
    return(genomic2plot(chr=chr, x=x, y=y, data.panel=data.panel, genome=genome, plot.params=plot.params))
  },
    ideogramMid=function(chr) {
      return(ideoMid(chr=chr, genome=genome, plot.params=plot.params))
  },
    chr.height=chrHeight(plot.params)    
  ))
}




####################################################################################
#
#    Plot Type 1 - Horizontal with Data Above the Ideogram  (BEGIN)
#
####################################################################################

getIdeogramMidY_1HorizDataAboveIdeogram <- function(chr, genome, plot.params) {
  pp <- plot.params
  chr.height <- getChrHeight_1HorizDataAboveIdeogram(pp)
  chr.names <- GenomeInfoDb::seqlevels(genome)
  chrs <- c(length(chr.names):1)
  names(chrs) <- chr.names
  chr.num <- chrs[chr]
  return(pp$bottommargin + (chr.num - 1) * chr.height + pp$ideogramheight/2)
}

getChrHeight_1HorizDataAboveIdeogram <- function(pp) {
  chr.height <- pp$ideogramheight + pp$data1inmargin + pp$data1height + pp$data1outmargin
  return(chr.height)
}

#Build the function mapping genomic regions into plotting coordinates
genomic2plot_1HorizDataAboveIdeogram <- function(chr=NULL, x=NULL, y=NULL, data.panel=1,
                                                 genome, plot.params) {
  
  if(is.null(data.panel)) data.panel <- 1
  if(data.panel != 1) {
    data.panel <- 1 #This plot type has only one data.panel
    warning("Specified data.panel does not exist. Plotting in data.panel 1")
  }
  
  pp <- plot.params
  
  genome.width <- 1 - pp$leftmargin - pp$rightmargin
  max.chr.len <- max(end(genome)) - min(start(genome))
    
  if(!is.null(x)) {
    x.plot <- pp$leftmargin + (x/max.chr.len)*genome.width
  } else{
    x.plot <- NULL
  }
  
  if(is.null(chr)) {
    y.plot <- NULL 
  } else {
    #chrs.y <- getChrLowestY(chr, genome, pp) 
    chrs.y <- getIdeogramMidY_1HorizDataAboveIdeogram(chr, genome, pp) 
    if(is.null(y)) { #if y is null, set y to the middle of the ideogram
      y.plot <- chrs.y
    } else { #Return the appropiate plot.y for the given original y
      datayrange <- pp$data1max - pp$data1min 
      yscaled <- ((y - pp$data1min) / datayrange) * pp$data1height
      y.plot <- chrs.y + pp$ideogramheight/2 + pp$data1inmargin + yscaled
    }
  }   
  return(list(x=x.plot, y=y.plot))
}

  

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#    Plot Type 1 - Horizontal with Data Above the Ideogram  (END)
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





####################################################################################
#
#    Plot Type 2 - Horizontal with Data Above and Below the Ideogram  (BEGIN)
#
####################################################################################

getIdeogramMidY_2HorizDataAboveAndBelowIdeogram <- function(chr, genome, plot.params) {
  pp <- plot.params
  chr.height <- getChrHeight_2HorizDataAboveAndBelowIdeogram(pp)
  chr.names <- GenomeInfoDb::seqlevels(genome)
  chrs <- c(length(chr.names):1)
  names(chrs) <- chr.names
  chr.num <- chrs[chr]
  return(pp$bottommargin + (chr.num - 1) * chr.height +
         pp$data2outmargin + pp$data2height + pp$data2inmargin + pp$ideogramheight/2)
}

getChrHeight_2HorizDataAboveAndBelowIdeogram <- function(pp) {
  chr.height <- (pp$data2outmargin + pp$data2height + pp$data2inmargin
                 + pp$ideogramheight
                 + pp$data1inmargin + pp$data1height + pp$data1outmargin)
  return(chr.height)
}

#Build the function mapping genomic regions into plotting coordinates
genomic2plot_2HorizDataAboveAndBelowIdeogram <- function(chr=NULL, x=NULL, y=NULL, data.panel,
                                                         genome, plot.params) {
  
  pp <- plot.params
  
  genome.width <- 1 - pp$leftmargin - pp$rightmargin
  max.chr.len <- max(end(genome)) - min(start(genome))
  
  if(!is.null(x)) {
    x.plot <- pp$leftmargin + (x/max.chr.len)*genome.width
  } else{
    x.plot <- NULL
  }
  
  if(is.null(chr)) {
    y.plot <- NULL 
  } else {
    #chrs.y <- getChrLowestY(chr, genome, pp) 
    chrs.y <- getIdeogramMidY_2HorizDataAboveAndBelowIdeogram(chr, genome, pp) 
    if(is.null(y)) { #if y is null, set y to the bottom of the chromosome
      y.plot <- chrs.y
    } else { #Return the appropiate plot.y for the given original y
      if(data.panel == 1) {
        datayrange <- pp$data1max - pp$data1min 
        yscaled <- ((y - pp$data1min) / datayrange) * pp$data1height
        y.plot <- chrs.y + pp$ideogramheight/2 + pp$data1inmargin + yscaled
      } else {
        if(data.panel == 2) {
          datayrange <- pp$data2max - pp$data2min 
          yscaled <- ((y - pp$data2min) / datayrange) * pp$data2height
          y.plot <- chrs.y - pp$ideogramheight/2 - pp$data2inmargin - yscaled
        } else {
          stop("Invalid data.panel")
        }
      }
    }
  }   
  return(list(x=x.plot, y=y.plot))
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#    Plot Type 2 - Horizontal with Data Above and Below the Ideogram  (END)
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





####################################################################################
#
#    Plot Type 4 - Vertical with Data Above and Below the Ideogram  (BEGIN)
#
####################################################################################

getIdeogramMidY_4VerticalDataAboveAndBelowIdeogram <- function(chr, genome, plot.params) {
  pp <- plot.params
  chr.height <- getChrHeight_4VerticalDataAboveAndBelowIdeogram(pp)
  chr.names <- GenomeInfoDb::seqlevels(genome)
  chrs <- c(length(chr.names):1)
  names(chrs) <- chr.names
  chr.num <- chrs[chr]
  return(pp$bottommargin + (chr.num - 1) * chr.height +
           pp$data2outmargin + pp$data2height + pp$data2inmargin + pp$ideogramheight/2)
}

getChrHeight_4VerticalDataAboveAndBelowIdeogram <- function(pp) {
  chr.height <- (pp$data2outmargin + pp$data2height + pp$data2inmargin
                 + pp$ideogramheight
                 + pp$data1inmargin + pp$data1height + pp$data1outmargin)
  return(chr.height)
}

#Build the function mapping genomic regions into plotting coordinates
genomic2plot_4VerticalDataAboveAndBelowIdeogram <- function(chr=NULL, x=NULL, y=NULL, 
                                                            data.panel, genome, plot.params) {
  
  pp <- plot.params
  
  genome.width <- 1 - pp$leftmargin - pp$rightmargin
  max.chr.len <- max(end(genome)) - min(start(genome))
  
  if(!is.null(x)) {
    x.plot <- pp$leftmargin + (x/max.chr.len)*genome.width
  } else{
    x.plot <- NULL
  }
  
  if(is.null(chr)) {
    y.plot <- NULL 
  } else {
    #chrs.y <- getChrLowestY(chr, genome, pp) 
    chrs.y <- getIdeogramMidY_4VerticalDataAboveAndBelowIdeogram(chr, genome, pp) 
    if(is.null(y)) { #if y is null, set y to the bottom of the chromosome
      y.plot <- chrs.y
    } else { #Return the appropiate plot.y for the given original y
      if(data.panel == 1) {
        datayrange <- pp$data1max - pp$data1min 
        yscaled <- ((y - pp$data1min) / datayrange) * pp$data1height
        y.plot <- chrs.y + pp$ideogramheight/2 + pp$data1inmargin + yscaled
      } else {
        if(data.panel == 2) {
          datayrange <- pp$data2max - pp$data2min 
          yscaled <- ((y - pp$data2min) / datayrange) * pp$data2height
          y.plot <- chrs.y - pp$ideogramheight/2 - pp$data2inmargin - yscaled
        } else {
          stop("Invalid data.panel")
        }
      }
    }
  }   
  return(list(x=y.plot, y=x.plot))
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#    Plot Type 4 - Vertical with Data Above and Below the Ideogram  (END)
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


