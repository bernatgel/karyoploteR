#' plotKaryotype
#' 
#' @description 
#' 
#' Create a new empty plot with a karyotype (the chromosome ideograms and chromosome names).
#' 
#' @details 
#'  
#'  This is the main function of \code{karyoploteR}. It creates the basic empty plot with 
#'  the chromosome ideograms and returns the karyoplot object needed for all other plotting 
#'  functions. Both the basic plotting parameters (margins, sizes, etc.) and the specific
#'  plotting functions for the ideograms and chromosome labels are customizable. 
#'  In particular, passing in a \code{plot.params} object specifies the basic plotting 
#'  parameters to use and the \code{ideogram.plotter} and \code{labels.plotter} parameters 
#'  can be used to specify custom plotting functions for the ideogram and the chromosome 
#'  labels. It is also possible to specify the genome and a list with the chromosomes to
#'  be plotted. 
#'  
#'  The \code{plot.type} parameter specifies the type of karyoplot to create: the number
#'  and positions of the data panels respect to the ideograms: 
#'  \itemize{
#'    \item \code{plot.type=1}  Horizontal ideograms with a single data panel above them
#'    \item \code{plot.type=2}  Horizontal ideograms with a two data panels, one above and one below them
#'  }
#'  
#'  More plot types are expected to come in the near future.
#'  
#' 
#' @usage plotKaryotype(genome="hg19", plot.type=1, ideogram.plotter=kpAddCytobands, labels.plotter=kpAddChromosomeNames, chromosomes="canonical", cytobands=NULL, plot.params=NULL, use.cache=TRUE, main=NULL, ...)
#' 
#' @param genome    The genome to plot. It can be either a UCSC style genome name (hg19, mm10, etc), a GRanges object with the chromosomes as ranges or in general any genome specification accepted by \code{\link[regioneR]{getGenomeAndMask}}. (defaults to "hg19")
#' @param plot.type    The orientation of the ideogram and placing of the data panels. Values explained above.. (defaults to 1)
#' @param ideogram.plotter    The function to be used to plot the ideograms. Only one function is included with the package, \code{kpAddCytobands}, but it is possible to create custom ones. If NULL, no ideograms are plotted. (defaults to \code{kpAddCytobands})
#' @param labels.plotter    The function to be used to plot the labels identifying the chromosomes. Only one function is included with the package, \code{kpAddChromosomeNames}, but it is possible to create custom ones. If NULL, no labels are plotted. (defaults to \code{kpAddChromosomeNames})
#' @param chromosomes    The chromosomes to plot. Can be either a vector of chromosome names or a chromosome group name ("canonical", "autosomal", "all"). (defaults to "canonical")
#' @param cytobands    A GRanges object specifying the positions and types of the cytobands. If NULL, the cytobands are recovered from the package cache or downloaded from UCSC. If empty, no cytobands will be plotted. (defaults to NULL)
#' @param plot.params    An object obtained from \code{\link{getDefaultPlotParams}} and possibly modified, containing the basic plotting parameters. If NULL, the defaults parameters will be used. (defaults to NULL)
#' @param use.cache    \code{karyoploteR} has a small cache with the chromosome names and lengths and the cytobands for a handful of organisms so it's not needed to retrieve them from databses or \code{BSGenomes} objects. Set this parameter to FALSE to ignore the cache. (defaults to TRUE, use the cache)
#' @param main    The text to be used as the title of the plot. NULL produces no title. (defaults to NULL)
#' @param ...    The ellipsis can be used to pass in any additional parameter accepted by the internal functions used.
#' 
#' @return
#' 
#' The \code{KaryoPlot} object needed by the plotting functions.
#' 
#' @seealso \code{\link{getDefaultPlotParams}}, \code{\link{kpPoints}}
#' 
#' @examples
#'  
#'  set.seed(1000)
#' 
#' rand.data <- createRandomRegions(genome="hg19", nregions=10000, length.mean=1, 
#'                                  length.sd=0, mask=NA, non.overlapping=TRUE)
#' mcols(rand.data) <- data.frame(y=rnorm(n=10000, mean = 0.5, sd=0.1))
#' 
#' #The simplest way, with all default parameters
#' kp <- plotKaryotype()
#' kpPoints(kp, rand.data, pch=".")
#' 
#' #Or we can plot only a few chromosomes, with 2 data panels
#' kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 2)
#' kpDataBackground(kp, data.panel = 1, color = "lightgreen")
#' kpDataBackground(kp, data.panel = 2, color = "lightblue")
#' kpPoints(kp, rand.data, pch=".", data.panel = 1)
#' kpPoints(kp, rand.data, pch=".", data.panel = 2)
#' 
#' #Or we can use a different organism, 
#' kp <- plotKaryotype(genome = "mm10")
#' kp <- plotKaryotype(genome = "dm6")
#' 
#' # Or we can change the plotting parameters. In this case, to create a smaller ideogram
#' # and smaller data panel below it
#' plot.params <- getDefaultPlotParams(plot.type=2)
#' plot.params$ideogramheight <- 5
#' plot.params$data2height <- 50
#' 
#' kp <- plotKaryotype(chromosomes = c("chr1", "chr2"), plot.type = 2, plot.params = plot.params)
#' kpDataBackground(kp, data.panel = 1, color = "lightgreen")
#' kpDataBackground(kp, data.panel = 2, color = "lightblue")
#' kpPoints(kp, rand.data, pch=".", data.panel = 1)
#' kpPoints(kp, rand.data, pch=".", data.panel = 2)
#' 
#' #Or we can remove the cytobands, passing an empty GRanges object
#' kp <- plotKaryotype(cytobands = GRanges())
#' 
#' #Or remove the chromosome labels
#' kp <- plotKaryotype(labels.plotter = NULL)
#' kpPoints(kp, rand.data, pch=".")
#' 
#' #In addition, it's possible to use maggrittr piping to chain the plotting calls
#' library(magrittr)
#' kp <- plotKaryotype() %>%
#'    kpDataBackground(color = "lightgreen") %>%
#'    kpPoints(rand.data, pch=".")
#' 
#' 
#' @import regioneR
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom S4Vectors runLength runValue
#' @importFrom memoise memoise
#' @importFrom rtracklayer ucscTableQuery
#' @importFrom biovizBase getBioColor
#' @importFrom biovizBase getIdeogram
#' @import methods
#' 
#' @export plotKaryotype
#' 


plotKaryotype <- function(genome="hg19", plot.type=1, ideogram.plotter=kpAddCytobands,
                          labels.plotter=kpAddChromosomeNames, chromosomes="canonical",
                          cytobands=NULL, plot.params=NULL, use.cache=TRUE, main=NULL, ...) {
  
  #check required parameters...
  
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(plot.type)
  }
  
  #Prepare the genome and filter the chromosomes as required
  #Get the genome
  #If the user has given us a valid GRanges as genome, use it directly
  #if it's something else, try with the cache or rely on regioneR::getGenomeAndMask
  if(is(genome, "GRanges")) {
    gr.genome <- genome
  } else {
    gr.genome <- NULL
    if(is(genome, "character") & use.cache==TRUE) { #Get the genome from the cache, if available
      if(genome %in% names(data.cache[["genomes"]])) {
        gr.genome <- data.cache[["genomes"]][[genome]]
      }
    }
    if(is.null(gr.genome)) {
      gr.genome <- getGenomeAndMask(genome=genome, mask=NA)$genome
    }
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
        message("There was an error when filtering the chromosomes. Using the unfiltered genome. \n", e)
      })     
    }
    # else Do not filter the chromosomes. If the genome is completely specified (not a character).
  }
  
  #Get the CytoBands if needed
  if(is.null(cytobands)) {
    if(is.character(genome)) {
      cytobands <- getCytobands(genome)
      #if there are cytobands, filter the cytobands using the current genome
      if(!is.null(cytobands) && length(cytobands)>0) {
        cytobands <- GenomeInfoDb::keepSeqlevels(cytobands, value=GenomeInfoDb::seqlevels(gr.genome), pruning.mode="coarse")
      }
    } else {
      message("No valid genome specified and no cytobands provided. No cytobands will be passed to the ideogram plotter.")
    }
  }

  
  
  #Create the KaryotypePlot Object that can be used to plot additional data onto the karyotype
    kp <- list()
    class(kp) <- "KaryoPlot"
    kp$plot.params <- plot.params
    if(is.character(genome)) {
      kp$genome.name <- genome
    } else {
      kp$genome.name <- "custom"
    }
    kp$chromosomes <- as.character(GenomeInfoDb::seqlevels(gr.genome))
    kp$chromosome.lengths <- setNames(end(gr.genome), seqnames(gr.genome))
    kp$genome <- gr.genome
    kp$cytobands <- cytobands
    kp$plot.type <- plot.type
    
    #Get the Coordinates Change Function to be used in this plot
    #coordChangeFunctions <- getCoordChangeFunctions(plot.type = plot.type, genome = gr.genome, plot.params = plot.params)
    coordChangeFunctions <- getCoordChangeFunctions(karyoplot=kp)
    kp$coord.change.function <- coordChangeFunctions$coordChangeFunction
    kp$ideogram.mid <- coordChangeFunctions$ideogramMid
    kp$chromosome.height <- coordChangeFunctions$chr.height
    
    
  #Remove all margins around the plot to take complete control of the available space
    kp$graphical.par <- list()
    kp$graphical.par$old.par <- graphics::par(no.readonly = TRUE)
    graphics::par(mar=c(0,0,0,0)+0.1)
    
    kp$beginKpPlot <- function() {
      graphics::par(kp$graphical.par$new.par)  
    }
    kp$endKpPlot <- function() {
      graphics::par(kp$graphical.par$old.par)  
    }
    
    on.exit(kp$endKpPlot())
  
  
  
  
  #Create the plot
  #TODO: Manage the specification of the y lab and xlab
    pp <- plot.params
    if(plot.type %in% c(1,2)) {
      xlim <- c(0, 1)
      ylim <- c(0, pp$bottommargin + length(gr.genome)*kp$chromosome.height + pp$topmargin)
    } else {
      if(plot.type %in% c(3,4,5)) {
        xlim <- c(0, 1)
        ylim <- c(0, pp$bottommargin + kp$chromosome.height + pp$topmargin)
      }
    }
  
    #create an empty plot
    graphics::plot(0, type="n", xlim=xlim, ylim=ylim, axes=FALSE, ylab="", xlab="", xaxs="i", yaxs="i")
  
  
    #Get the limits of the plot from the graphical device
    kp$plot <- list()
      p <- graphics::par("usr")
      kp$plot$xmin <- p[1]
      kp$plot$xmax <- p[2]
      kp$plot$ymin <- p[3]  
      kp$plot$ymax <- p[4]
 
  kp$graphical.par$new.par <- graphics::par(no.readonly = TRUE) #Remember the parameters used

  #Finally, plot the ideograms and labels, if needed
  #And plot the ideogram
  if(!is.null(ideogram.plotter)) {
    ideogram.plotter(kp, ...)
  }    

  #Plot the Chromosome Labels
    if(!is.null(labels.plotter)) {
      labels.plotter(kp, ...)
    }  
  
  #Add the main title
  if(!is.null(main)) {
    kpAddMainTitle(kp, main, ...)
  }
  
  
 
  return(kp)
}
