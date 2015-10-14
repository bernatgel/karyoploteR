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



library(regioneR)
library(rtracklayer)
library(memoise)


plotKaryotype <- function(genome="hg19", ideogram.plotter=plotCytobands, labels.plotter=plotChromosomeNames, chromosomes="canonical", cytobands=NULL, plot.params=NULL, ...) {
  
  #check required parameters...
  
  if(is.null(plot.params)) {
    #TODO: Move to a function that returns the default parameters
    plot.params <- list(xleftmargin=0.1, xrightmargin=0.05, ytopmargin=100, ybottommargin=100,
                   yabovemargin=5, ybelowmargin=5, ydataheight=200, ideogramheight=50,
                   dataymin=0, dataymax=100)
  }
  
  #Prepare the genome and filter the chromosomes as required
  gr.genome <- getGenomeAndMask(genome=genome, mask=NA)$genome
  if(!is.null(chromosomes) & chromosomes != "all") {
    if(is.character(genome)) {
      if(chromosomes %in% c("canonical", "autosomal")) {
        gr.genome <- filterChromosomes(gr.genome, organism=genome, chr.type=chromosomes)
      } else {
        gr.genome <- filterChromosomes(gr.genome, keep.chr=chromosomes)
      }
    } else {
      #Do not filter the chromosomes. If the genome is completely specified.
    } 
  }
  print(gr.genome)
  #Get the CytoBands if needed
  if(is.null(cytobands)) {
    if(is.character(genome)) {
      cytobands <- getCytobands(genome)
      #Filter the cytobands using the current genome
      cytobands <- keepSeqlevels(cytobands, value=seqlevels(gr.genome))
      
    } else {
      message("No valid genome specified and no cytobands provided. No cytobands will be passed to the ideogram plotter.")
    }
  }

  #Get the Coordinates Change Function to be used in this plot
  coordChangeFunction <- getCoordChangeFunction(gr.genome, plot.params)
  
  
  #Create the KaryotypePlot Object that can be used to plot additional data onto the karyotype
    kp <- list()
    class(kp) <- "KaryoPlot"
    kp$plot.params <- plot.params
    kp$coord.change.function <- coordChangeFunction
    if(is.character(genome)) {
      kp$genome.name <- genome
    } else {
      kp$genome.name <- "custom"
    }
    kp$chromosomes <- as.character(seqlevels(gr.genome))
    kp$genome <- gr.genome
    kp$cytobands <- cytobands
  
  
  #Create the plot
  #TODO: Manage the specification of the y lab and xlab
    pp <- plot.params
    xlim <- c(0, 1)
    chr.height <- pp$ybelowmargin + pp$ideogramheight + pp$yabovemargin + pp$ydataheight
    ylim <- c(0, pp$ybottommargin + length(gr.genome)*chr.height + pp$ytopmargin)
    
    #create an empty plot
    #TODO: Should we remove any margin around the plot?
    par(mar=c(0,0,0,0)+0.1)
    plot(0, type="n", xlim=xlim, ylim=ylim, axes=FALSE, ylab="", xlab="", xaxs="i", yaxs="i")
  
    #Get the limits of the plot from the graphical device
    kp$plot <- list()
      p <- par("usr")
      kp$plot$xmin <- p[1]
      kp$plot$xmax <- p[2]
      kp$plot$ymin <- p[3]  
      kp$plot$ymax <- p[4]  
  
  #Plot the Chromosome Labels
    if(!is.null(labels.plotter)) {
      labels.plotter(kp, ...)
    }  
 
  #And plot the ideogram
  if(!is.null(ideogram.plotter)) {
    ideogram.plotter(kp, ...)
  }  
  
  return(kp)
}

genome <- data.frame(c("A", "B"), c(0, 0), c(20000000, 15000000))
kp <- plotKaryotype(genome="hg19")
kp$plot
kp$coord.change.function(x=0)
  
available.genomes()
  
text(x=coordChangeFunction(chr="chr17", x=29000000)$x, y=coordChangeFunction(chr="chr17", y=100)$y, labels="NF1")
  

kpPlotRegions(kp, createRandomRegions(nregions=20, length.mean=10000000, length.sd=10000000, mask=NA)) 
kpPlotRegions(kp, createRandomRegions(nregions=20, length.mean=10000000, length.sd=10000000, mask=NA), y=10) 
kpPlotRegions(kp, createRandomRegions(nregions=20, length.mean=10000000, length.sd=10000000, mask=NA), y=70)   

abline(v=kp$coord.change.function(x=max(end(kp$genome)))$x)
abline(v=kp$coord.change.function(x=min(start(kp$genome)))$x)
abline(v=1)
abline(v=0)

kpPlotRegions(coordChangeFunction, regions)

kp$plot


