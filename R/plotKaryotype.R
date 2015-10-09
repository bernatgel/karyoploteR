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


plotKaryotype <- function(genome="hg19", ideogram.plotter, labels.plotter=plotChromosomeNames, chromosomes="canonical", cytobands=NULL, plot.params=NULL, ...) {
  
  #check required parameters...
  
  if(is.null(plot.params)) {
    #TODO: Move to a function that returns the default parameters
    plot.params <- list(xleftmargin=100, xrightmargin=20, ytopmargin=20, ybottommargin=20,
                   yabovemargin=5, ybelowmargin=5, ydataheight=100, ideogramheight=50,
                   dataymin=0, dataymax=100)
  }
  
  #Prepare the genome and filter the chromosomes as required
  gr.genome <- getGenomeAndMask(genome=genome, mask=NA)$genome
  if(!is.null(chromosomes) & chromosomes != "all") {
    if(is.character(genome) & chromosomes %in% c("canonical", "autosomal")) {
      gr.genome <- filterChromosomes(gr.genome, organism=genome, chr.type=chromosomes)
    } else {
      gr.genome <- filterChromosomes(gr.genome, keep.chr=chromosomes)
    }
  }
  
  #Get the CytoBands if needed
  if(is.null(cytobands)) {
    if(is.character(genome)) {
      cytobands <- getCytobands(genome)
    } else {
      message("No valid genome specified and no cytobands provided. No cytobands will be passed to the ideogram plotter.")
    }
  }

  #Get the Coordinates Change Function to be used in this plot
  coordChangeFunction <- getCoordChangeFunction(gr.genome, plot.params)
  
  #Create the plot
  #TODO: Manage the specification of the y lab and xlab
    pp <- plot.params
    xlim <- c(0, pp$xleftmargin + coordChangeFunction(x=max(end(gr.genome)))$x + pp$xrightmargin)
    chr.height <- pp$ybelowmargin + pp$ideogramheight + pp$yabovemargin + pp$ydataheight
    ylim <- c(0, pp$ybottommargin + length(gr.genome)*chr.height + pp$ytopmargin)
    
    #create an empty plot
    #TODO: Should we remove any margin around the plot?
    par(mar=c(0,0,0,0)+0.1)
    plot(0, type="n", xlim=xlim, ylim=ylim, axes=FALSE, ylab="", xlab="")
  
    
  #Plot the Chromosome Labels
    if(!is.null(labels.plotter)) {
      labels.plotter(coordChangeFunction, gr.genome, plot.params, ...)
    }  
  
  #Plot the ideograms
  
  cytobands <- keepSeqlevels(cytobands, value=seqlevels(gr.genome))
  
  ybottom <- coordChangeFunction(as.character(seqnames(cytobands)))$y + pp$ybelowmargin
  ytop <- coordChangeFunction(as.character(seqnames(cytobands)))$y + pp$ybelowmargin + pp$ideogramheight
    
  xleft <- coordChangeFunction(x=start(cytobands))$x
  xright <- coordChangeFunction(x=end(cytobands))$x

  col <- do.call(c, colorTable[cytobands$gieStain])
  
  rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)
}



colorTable <- list(gneg="gray90",
                 gpos25="gray75",
                 gpos33="gray66",
                 gpos50="gray50",
                 gpos66="gray33",
                 gpos75="gray25",
                 pos100=
                 stalk="blue", #repetitive areas
                 acen="darkred", #centromeres
                 gvar="black")

  
colors = {
  'gpos100' : (0/255.0,0/255.0,0/255.0),
  'gpos' : (0/255.0,0/255.0,0/255.0),
  'gpos75' : (130/255.0,130/255.0,130/255.0),
  'gpos66' : (160/255.0,160/255.0,160/255.0),
  'gpos50' : (200/255.0,200/255.0,200/255.0),
  'gpos33' : (210/255.0,210/255.0,210/255.0),
  'gpos25' : (200/255.0,200/255.0,200/255.0),
  'gvar' : (220/255.0,220/255.0,220/255.0),
  'gneg' : (255/255.0,255/255.0,255/255.0),
  'acen' : (217/255.0,47/255.0,39/255.0),
  'stalk' : (100/255.0,127/255.0,164/255.0),
}

  
  text(x=coordChangeFunction(chr="chr17", x=29000000)$x, y=coordChangeFunction(chr="chr17", y=100)$y, labels="NF1")
  
  kpPlotRegions <- function(coordChangeFunction, regions, y=50, height=20, ...) {
    xleft <- coordChangeFunction(x=start(regions))$x
    xright <- coordChangeFunction(x=end(regions))$x
    ytop <- coordChangeFunction(chr=as.character(seqnames(regions)), y=rep(y+height/2, length(regions)))$y
    ybottom <- coordChangeFunction(chr=as.character(seqnames(regions)), y=rep(y-height/2, length(regions)))$y
    rect(xleft=xleft, xright=xright, ytop=ytop, ybottom=ybottom, col="red")
  }
  
  



