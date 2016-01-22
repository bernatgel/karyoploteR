#' plotKaryotype
#' 
#' @description 
#' 
#' Create a new plot with a karyotype and nothing more and return the coordinates change function
#' 
#' @details 
#'  
#' 
#' @usage plotKaryotype(genome="hg19", ...)
#' 
#' @param genoms
#' 
#' @return
#' The coordinates change function to be used to plot additional data onto the karyotype 
#' 
#' 
#' @seealso \code{\link{getMask}}, \code{\link{getGenomeAndMask}}, \code{\link{characterToBSGenome}}, \code{\link{maskFromBSGenome}}, \code{\link{emptyCacheRegioneR}}
#' 
#' @examples
#'  
#' @export plotKaryotype
#' 


plotKaryotype <- function(genome="hg19", plot.type=1, ideogram.plotter=plotCytobands, labels.plotter=plotChromosomeNames, chromosomes="canonical", cytobands=NULL, plot.params=NULL, use.cache=TRUE, ...) {
  
  #check required parameters...
  
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(plot.type)
  }
  
  #Prepare the genome and filter the chromosomes as required
  #Get the genome
    gr.genome <- NULL
    if(use.cache) { #Get the genome from the cache, if available
      if(genome %in% names(data.cache[["genomes"]])) {
        gr.genome <- data.cache[["genomes"]][[genome]]
      }
    }
    if(is.null(gr.genome)) {
      gr.genome <- getGenomeAndMask(genome=genome, mask=NA)$genome
    }
  #And filter it
  if(!is.null(chromosomes) && chromosomes != "all") {
    if(is.character(genome)) {
      tryCatch(expr={
        if(length(chromosomes)==1 && (chromosomes %in% c("canonical", "autosomal"))) {
          gr.genome <- filterChromosomes(gr.genome, organism=genome, chr.type=chromosomes)
        } else {
          gr.genome <- filterChromosomes(gr.genome, keep.chr=chromosomes)
        }
      }, error=function(e) {
        message(paste0("There was an error when filtering the chromosomes. Using the unfiltered genome. \n", e))
      })     
    }
    # else Do not filter the chromosomes. If the genome is completely specified (not a character).
  }
  
  #Get the CytoBands if needed
  if(is.null(cytobands)) {
    if(is.character(genome)) {
      cytobands <- getCytobands(genome)
      #if there are sytobands, filter the cytobands using the current genome
      if(!is.null(cytobands) && length(cytobands)>0) {
        cytobands <- keepSeqlevels(cytobands, value=seqlevels(gr.genome))
      }
    } else {
      message("No valid genome specified and no cytobands provided. No cytobands will be passed to the ideogram plotter.")
    }
  }

  #Get the Coordinates Change Function to be used in this plot
  coordChangeFunctions <- getCoordChangeFunctions(plot.type, gr.genome, plot.params)
  
  
  #Create the KaryotypePlot Object that can be used to plot additional data onto the karyotype
    kp <- list()
    class(kp) <- "KaryoPlot"
    kp$plot.params <- plot.params
    kp$coord.change.function <- coordChangeFunctions$coorChangeFunction
    kp$ideogram.mid <- coordChangeFunctions$ideogramMid
    kp$chromosome.height <- coordChangeFunctions$chr.height
    if(is.character(genome)) {
      kp$genome.name <- genome
    } else {
      kp$genome.name <- "custom"
    }
    kp$chromosomes <- as.character(seqlevels(gr.genome))
    kp$genome <- gr.genome
    kp$cytobands <- cytobands
  
  
  #Remove all margins around the plot to take complete control of the available space
    kp$graphical.par <- list()
    kp$graphical.par$old.par <- par(no.readonly = TRUE)
    par(mar=c(0,0,0,0)+0.1)
    
    kp$beginKpPlot <- function() {
      par(kp$graphical.par$new.par)  
    }
    kp$endKpPlot <- function() {
      par(kp$graphical.par$old.par)  
    }
    
    on.exit(kp$endKpPlot())
  
  
  
  
  #Create the plot
  #TODO: Manage the specification of the y lab and xlab
    pp <- plot.params
    if(plot.type %in% c(1,2,3)) {
      xlim <- c(0, 1)
      ylim <- c(0, pp$bottommargin + length(gr.genome)*kp$chromosome.height + pp$topmargin)
    } else {
      ylim <- c(0, 1)
      xlim <- c(0, pp$bottommargin + length(gr.genome)*kp$chromosome.height + pp$topmargin)
    }
  
    #create an empty plot
    plot(0, type="n", xlim=xlim, ylim=ylim, axes=FALSE, ylab="", xlab="", xaxs="i", yaxs="i")
  
  
    #Get the limits of the plot from the graphical device
    kp$plot <- list()
      p <- par("usr")
      kp$plot$xmin <- p[1]
      kp$plot$xmax <- p[2]
      kp$plot$ymin <- p[3]  
      kp$plot$ymax <- p[4]
  
  #And plot the ideogram
  if(!is.null(ideogram.plotter)) {
    ideogram.plotter(kp, ...)
  }    
  
  #Plot the Chromosome Labels
    if(!is.null(labels.plotter)) {
      labels.plotter(kp, ...)
    }  
 
  
  kp$graphical.par$new.par <- par(no.readonly = TRUE) #Remember the parameters used
  
  return(kp)
}
