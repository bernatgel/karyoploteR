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
#'    \item \code{plot.type=2}  Horizontal ideograms with two data panels, one above and one below them
#'    \item \code{plot.type=3}  Horizontal ideograms with all chromosomes in a single line with two data panels, one above and one below them
#'    \item \code{plot.type=4}  Horizontal ideograms with all chromosomes in a single line with one data panel above
#'    \item \code{plot.type=5}  Horizontal ideograms with all chromosomes in a single line with one data panel below them
#'    \item \code{plot.type=6}  Horizontal ideograms with NO data panels. Only plotting in the ideograms is possible.
#'    \item \code{plot.type=7}  Horizontal ideograms with all chromosomes in a single line with NO data panels. Only plotting in the ideograms is possible.
#'  }
#'  
#' There's more information at the \url{https://bernatgel.github.io/karyoploter_tutorial/}{karyoploteR tutorial}.
#'  
#' 
#' @usage plotKaryotype(genome="hg19", plot.type=1, ideogram.plotter=kpAddCytobands, labels.plotter=kpAddChromosomeNames, chromosomes="auto", zoom=NULL, cytobands=NULL, plot.params=NULL, use.cache=TRUE, main=NULL, ...)
#' 
#' @param genome    The genome to plot. It can be either a UCSC style genome name (hg19, mm10, etc), a BSgenome, a Seqinfo object, a GRanges object with the chromosomes as ranges or in general any genome specification accepted by \code{\link[regioneR]{getGenomeAndMask}}. (defaults to "hg19")
#' @param plot.type    The orientation of the ideogram and placing of the data panels. Values explained above.. (defaults to 1)
#' @param ideogram.plotter    The function to be used to plot the ideograms. Only one function is included with the package, \code{kpAddCytobands}, but it is possible to create custom ones. If NULL, no ideograms are plotted. (defaults to \code{kpAddCytobands})
#' @param labels.plotter    The function to be used to plot the labels identifying the chromosomes. Only one function is included with the package, \code{kpAddChromosomeNames}, but it is possible to create custom ones. If NULL, no labels are plotted. (defaults to \code{kpAddChromosomeNames})
#' @param chromosomes    The chromosomes to plot. Can be either a vector of chromosome names or a chromosome group name ("canonical", "autosomal", "all"). Setting it yo "auto" will select canonical for named genomes and no filtering for custom genomes. (defaults to "auto")
#' @param zoom   A GRanges object specifiyng a single region to zoom in or any format accepted by \code{regioneR::toGRanges}. If not NULL, it takes precedence over \code{chromosome} and only the zoomed in region is represented. If more than one region is present in the GRanges, only the first one is used. (defaults to NULL, do not zoom in and show the whole plot as specified by \code{genome} and \code{chromosomes})
#' @param cytobands    A GRanges object (or anything accepted by \code{\link[regioneR]{toGRanges}} function: bed file, data.frame...) specifying the positions and types of the cytobands. The type of the cytobands MUST be in a column named "gieStain" (as used by UCSC) with values such as 'gneg', 'gpos50', 'stalk', 'acen'... If NULL, the cytobands are recovered from the package cache or downloaded from UCSC if possible (it's not possible for custom genomes). If empty, no cytobands will be plotted. (defaults to NULL)
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
#' @importFrom IRanges overlapsAny
#' @import methods
#' @importFrom stats setNames
#' 
#' @export plotKaryotype
#' 


plotKaryotype <- function(genome="hg19", plot.type=1, ideogram.plotter=kpAddCytobands,
                          labels.plotter=kpAddChromosomeNames, chromosomes="auto",
                          zoom=NULL, cytobands=NULL, plot.params=NULL,
                          use.cache=TRUE, main=NULL, ...) {
  
  #Parameters Check
  #TODO: Finish checks
  #genome
  if(is.null(genome)) stop("genome cannot be NULL.")
  #plot.type
  
  
  #cytobands
  if(!is.null(cytobands)) {
    cytobands <- tryCatch(toGRanges(cytobands), error=function(e) {})
    if(!methods::is(cytobands, "GRanges")) stop("'cytobands' must be NULL, a GRanges or something accepted by regioneR::toGRanges")
  }
  
  #zoom
  if(!is.null(zoom)) {
    zoom <- regioneR::toGRanges(zoom)
    if(!methods::is(zoom, "GRanges")) stop("'zoom' must be NULL or a GRanges object")
    if(length(zoom)>1) {
      warning("The zoom parameter has more than one region. Only the first one will be used.")
      zoom <- zoom[1]
    }
  }
  
  #End parameters check
  
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(plot.type)
  }
  
  
  #Prepare the genome and filter the chromosomes as required
  #Get the genome
  #If the user has given us a valid GRanges as genome, use it directly
  #if it's something else, try with the cache or rely on regioneR::getGenomeAndMask
  gr.genome <- NULL
  genome.name <- NULL
  if(methods::is(genome, "GRanges")) {
    gr.genome <- genome
  } else {
    if(methods::is(genome, "BSgenome")) {
      genome <- seqinfo(genome)
    }
    if(methods::is(genome, "Seqinfo")) {
      gr.genome <- as(genome, "GRanges")
      genome.name <- genome(genome)[1]
    } else {
      if(methods::is(genome, "character") & use.cache==TRUE) { #Get the genome from the cache, if available
        if(genome %in% names(data.cache[["genomes"]])) {
          gr.genome <- data.cache[["genomes"]][[genome]]
          genome.name <- genome
        }
      }
    }
  }
  if(is.null(gr.genome)) {
    gr.genome <- tryCatch(regioneR::getGenomeAndMask(genome=genome, mask=NA)$genome,
                          error=function(e) {stop("It was not possible to identify or load the requested genome. ", e)}
    )
  }
  
  
  #Check the genome has no problems (repeated chromosomes, etc...)
  chr.names <- as.character(GenomeInfoDb::seqnames(gr.genome))
  if(any(duplicated(chr.names))) {
    stop(paste0("There are duplicate chromosome names in the genome. Chromosome names must be unique. Chromosome names are: ", paste0(chr.names, collapse = ", ")))
  }
  
  #Use the chromosome names as the names of the genome GRanges
  names(gr.genome) <- chr.names
  

  #Chromosome Filtering
  
  #If zoom is set, change the chromosome parameter to the name of the chomosome
  #we are zooming in, so everything else is automatically filtered out.
  if(!is.null(zoom)) {
    if(!IRanges::overlapsAny(zoom, gr.genome)) {
      stop("You are trying to set the zoom to a region not part of the current genome.")
    } else {
      chromosomes <- as.character(GenomeInfoDb::seqnames(zoom))
    }
  }
  
  #And filter it
  if(!is.null(chromosomes) && any(chromosomes != "all")) {
    
    if(length(chromosomes)==1 && (chromosomes %in% c("canonical", "autosomal", "auto"))) {
      if(!is.null(genome.name) && is.character(genome.name)) {
        if(chromosomes=="auto") chromosomes <- "canonical"   #Set it to canonical to perform the actual filtering
        tryCatch(expr={gr.genome <- filterChromosomes(gr.genome, organism=genome.name, chr.type=chromosomes)},
                 error=function(e) {
                   message("WARNING: There was an error when filtering the chromosomes and selecting only ", chromosomes, " chromosomes.  Falling back to using the unfiltered genome. \n", e)
                 })
      } else {
        if(chromosomes != "auto") { #If set to 'auto', say nothing about the filtering. The user has actually not requested any filtering on their GRanges.
          message("NOTE: It is only possible to filter chromosomes using named ",
                  "chromosome classes (i.e. 'canonical', 'autosomal') when the genome ",
                  "is specified by name (i.e. 'hg19', 'mm10'). Please, either ",
                  "use a genome specified by name or explicitly select the ",
                  "chromosomes to plot (i.e. chromosomes=c('chr1', 'chr2') ). ",
                  " Falling back to using the unfiltered genome.")
        }
      }
    } else {
      if(!all(chromosomes %in% as.character(seqnames(gr.genome)))) {
        message("NOTE: Not all requested chromosomes are part of the genome. Trying to filter as requested. ",
                "   * Requested Chromosomes: ", paste0(chromosomes, collapse = ", "),
                "   * Chromosomes in the genome: ", paste0(as.character(seqnames(gr.genome)), collapse = ", ")
        )
      }
      tryCatch(expr={gr.genome <- filterChromosomes(gr.genome, keep.chr=chromosomes)[chromosomes]}, #The selection by "[chromsomes]" ensures the order is the one specified and not the canonical in the genome
               error=function(e) {
                 message("WARNING: There was an error when filtering the chromosomes. Falling back to using the unfiltered genome. \n", e)
               }
      )
    }
  }
  
  #Check we still have a genome! Explanation: If the filtering step partially fails, we might end up with an empty genome
  if(length(gr.genome)==0) {
    stop("The genome has no chromosomes left after filtering. Cannot plot with no chromosomes.")
  }  

  
  #Get the CytoBands if needed
  if(is.null(cytobands)) {
    if(!is.null(genome.name) && all(is.character(genome.name))) {
      cytobands <- getCytobands(genome.name)
    } else {
      if(all(is.character(genome))) {
        cytobands <- getCytobands(genome)
      } 
    }
    #if there are cytobands, filter the cytobands using the current genome
    if(!is.null(cytobands) && length(cytobands)>0) {
      cytobands <- GenomeInfoDb::keepSeqlevels(cytobands, value=GenomeInfoDb::seqlevels(gr.genome), pruning.mode="coarse")
    }
  }

  #Create the KaryoPlot Object that can be used to plot additional data onto the karyotype
    kp <- list()
    class(kp) <- "KaryoPlot"
    kp$plot.params <- plot.params
    if(!is.null(genome.name)) {
      kp$genome.name <- genome.name
    } else {
      kp$genome.name <- "custom"
    }
    kp$chromosomes <- as.character(GenomeInfoDb::seqlevels(gr.genome))
    kp$chromosome.lengths <- stats::setNames(end(gr.genome), seqnames(gr.genome))
    kp$genome <- gr.genome
    kp$cytobands <- cytobands
    kp$plot.type <- plot.type
    
    
    #If zoom is NULL, set the plot.region to the whole genome. If it's not null, set it to the zoom region
    if(is.null(zoom)) {
      kp$plot.region <- kp$genome
      kp$zoom <- FALSE
    } else {
      kp$plot.region <- zoom
      names(kp$plot.region) <- as.character(seqnames(kp$plot.region))
      kp$zoom <- TRUE
    }
    
    
    
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
    pp <- plot.params
    if(plot.type %in% c(1,2,6)) {
      xlim <- c(0, 1)
      ylim <- c(0, pp$bottommargin + length(gr.genome)*kp$chromosome.height + pp$topmargin)
    } else {
      if(plot.type %in% c(3,4,5,7)) {
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
