#' kpPlotGenes
#' 
#' @description 
#' 
#' Plots rectangles along the genome representing the regions (or intervals) specified by a \code{GRanges} object
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of karyoploteR. It takes a \code{GRanges} object and
#'  plots its content. Overlapping regions can be stacked and the number of layers for overlapping regions can be set.
#'  In contrast with the low-level functions such as \code{\link{kpRect}}, it is not possible to specify the data using 
#'  independent numeric vectors and the function only takes in \code{GRanges}.
#'
#' @usage kpPlotGenes(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="black", border=NULL, avoid.overlapping=TRUE, num.layers=NULL, layer.margin=0.05, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the regions to plot.
# #removed as requested by the package reviewer. It can be any of the formats accepted by the \code{\link[regioneR]{toGRanges}} function from the package \href{http://bioconductor.org/packages/release/bioc/html/regioneR.html}{regioneR}.
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param col    (color) The background color of the regions. (defaults to black)
#' @param border    (color) The color used to draw the border of the regions. If NULL, no border is drawn. (defaults to NULL)
#' @param avoid.overlapping    (boolean) Whether overlapping regions should be drawn as stacks (TRUE) on drawing one occluding the other in a single layer (FALSE). (defaults to TRUE)
#' @param num.layers    (numeric) The number of layers the plotting space should be divided into to allow for plotting overlapping regions. The lotting region will be divided into this many pieces regardless if any overlapping regions actually exist. If NULL, the maximum number of regions overlapping a single point in the genome. (defaults to NULL)
#' @param layer.margin    (numeric) The blank space left between layers of regions. (defaults to 0.05)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#'  
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpSegments}}
#' 
#' @examples
#'  
#'  
#'  set.seed(1000)
#'  
#'  #Example 1: create 20 sets of non-overlapping random regions and plot them all. Add a coverage plot on top.
#'  kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#'  
#'  all.regs <- GRanges()
#'  
#'  nreps <- 20
#'  for(i in 1:nreps) {
#'    regs <- createRandomRegions(nregions = 100, length.mean = 10000000, length.sd = 1000000,
#'                                non.overlapping = TRUE, genome = "hg19", mask=NA)
#'    all.regs <- c(all.regs, regs)
#'    kpPlotGenes(kp, regs, r0 = (i-1)*(0.8/nreps), r1 = (i)*(0.8/nreps), col="#AAAAAA")
#'  }
#'  
#'  kpPlotCoverage(kp, all.regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
#'  kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
#'  
#'  
#'  #Example 2: Do the same with a single bigger set of possibly overlapping regions
#'  
#'  kp <- plotKaryotype("hg19", plot.type=1, chromosomes=c("chr1", "chr2"))
#'  
#'  regs <- createRandomRegions(nregions = 1000, length.mean = 10000000, length.sd = 1000000,
#'                              non.overlapping = FALSE, genome = "hg19", mask=NA)
#'                              
#'  kpPlotGenes(kp, regs, r0 = 0, r1 = 0.8, col="#AAAAAA")
#'  
#'  kpPlotCoverage(kp, regs, ymax = 20, r0=0.8,  r1=1, col="#CCCCFF")
#'  kpAxis(kp, ymin = 0, ymax= 20, numticks = 2, r0 = 0.8, r1=1)
#'  
#'  
#'  
#'@export kpPlotGenes


kpPlotGenes <- function(karyoplot, data, plot.transcripts=TRUE, plot.transcripts.structure=TRUE,
                        gene.margin=0.3, gene.col="black",
                        add.gene.name=TRUE, gene.names=NULL, gene.name.position="top", gene.name.cex=1, gene.name.col="black",
                        transcript.margin=0.5, transcript.col="black", 
                        add.transcript.name=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=0.6, transcript.name.col="black",
                        data.panel=1, r0=NULL, r1=NULL, col="black", 
                          border=NULL, avoid.overlapping=TRUE, num.layers=NULL,
                          layer.margin=0.05, clipping=TRUE, ...) {
  
  #karyoplot
    if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(missing(data)) stop("The parameter 'data' is required")
    if(!methods::is(data, "TxDb")) {
      if(!isValidData(data)) {
        stop("'data' must be either a TxDb object or a list with the required slots. See ?kpPlotGenes for more information")  
      }
    } 
  
  #If data is a TxDb object, build a list from it with the expected format
  if(methods::is(data, "TxDb")) {
    data <- tryCatch(makeDataFromTxDb(karyoplot, data, plot.transcripts, plot.transcripts.structure),
                     error=function(e) {"Error: There was an error extracting the information from the TxDb object."})
  }
  
  #if there's nothing to plot, return
  if(length(data$genes)==0) {
    #add an empty "latest.plot" information
    #TODO
    #and return
    invisible(karyoplot)
  }
  
  if(length(data$genes)>20 & plot.transcripts.structure==TRUE) {
    message("NOTE: Plotting many genes with detailed transcript structure may take a long time. You can set 'plot.transcripts' and 'plot.transcripts.structure' to FALSE to speed up the process reducing the detail level.")
  }
  
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  total.height <- 1
  
  if(plot.transcripts==FALSE) {
    #Plot only the genes, with one rectangle per gene
    genes.for.coverage <- data$genes
    strand(genes.for.coverage) <- "*"
    bins <- disjointBins(genes.for.coverage)
    num.layers <- max(bins)
    layer.height <- total.height/num.layers
    gene.height <- layer.height/(1+gene.margin)
    y0 <- layer.height * (bins-1)
    y1 <- y0 + gene.height
    kpRect(karyoplot, data=data$genes, y0=y0, y1=y1, col=gene.col)
    if(add.gene.name==TRUE) {
      gene.labels <- ifelse(!is.null(gene.names), gene.names[names(data$genes)], names(data$genes))
      plotNames(karyoplot, data=data$genes, y0=y0, y1=y1, labels=gene.labels, position=gene.name.position, col=gene.name.col, cex=gene.name.cex, r0=r0, r1=r1, data.panel=data.panel)
    }
    #TODO: Store the position of the genes in the latest.plot slot
  } else { #if plot.transcripts==TRUE
    #Compute the position of each transcript (and the height dedicated to each gene depending on the numbre of transcripts)
  
    #Treat all genes as if they have the same length as the gene, so each genes has its own "rectangular space" preserved
    num.transcripts.per.gene <- Map(length, as.list(data$transcripts))
    genes.for.coverage <- Reduce(c, Map(rep, as.list(data$genes), num.transcripts.per.gene))
    strand(genes.for.coverage) <- "*"
    max.coverage <- max(max(coverage(genes.for.coverage)))
    
    bin.height <- total.height/max.coverage
    transcript.height <- bin.height/(1+transcript.margin)
    
    #Once done, plot every transcript separately 
    #and position them using the binning functionality of Bioconductor
    bins <- disjointBins(genes.for.coverage)
    names(bins) <- names(genes.for.coverage)
    
    col.num <- 1
    last.gene <- ""
    last.gene.y0 <- 0
    last.gene.y1 <- 0
    
    transcript.in.gene <- 1
    for(nt in seq_len(length(bins))) {
      bin <- bins[nt]
      message(bin)
      if(last.gene != names(bin)) {
        #Plot the gene label if needed
        if(last.gene.y0 != last.gene.y1) {
          gene.labels <- ifelse(!is.null(gene.names), gene.names[last.gene], last.gene)
          plotNames(karyoplot, data=data$genes[last.gene], y0=last.gene.y0, y1=last.gene.y1, labels=gene.labels, position=gene.name.position, col=gene.name.col, cex=gene.name.cex, r0=r0, r1=r1, data.panel=data.panel)
          #TODO: Store the y0, y1 position of the gene into the latest plot object
        }
        last.gene.y0 <- last.gene.y1
        #update the counters and flags
        transcript.in.gene <- 1
        last.gene <- names(bin)
        col.num <- col.num+1
        tcol <- cols[col.num]
        message("new col: ", tcol)
        
      } else {
        transcript.in.gene <- transcript.in.gene + 1
      }
      
      transcript <- data$transcripts[[names(bin)]][transcript.in.gene]
      
      transcript.y0 <- bin.height*(bin-1)
      transcript.y1 <- transcript.y0 + transcript.height
      last.gene.y1 <- transcript.y1 #store the latest transcript y1 to later plot the gene.name if needed
      
      if(plot.transcripts.structure==FALSE) {
        kpRect(karyoplot, data=transcript, y0=transcript.y0, y1=transcript.y1, col=tcol, r0=r0, r1=r1, data.panel=data.panel)
        
        if(add.transcript.name==TRUE) {
          if(!is.null(transcript.names)) {
            transcript.labels <- transcript.names[transcript$tx_id]
          } else {
            transcript.labels <- transcript$tx_id
          }
          message(transcript.y0)
          plotNames(karyoplot, data=transcript, y0=transcript.y0, y1=transcript.y1, labels=transcript.labels, position=transcript.name.position, col=transcript.name.col, cex=transcript.name.cex, r0=r0, r1=r1, data.panel=data.panel)
        }
        
      } else { #if plot.transcript.structure==TRUE
        coding.exons <- data$coding.exons[[as.character(transcript$tx_id)]]
        non.coding.exons <- data$non.coding.exons[[as.character(transcript$tx_id)]]
        plotTranscript(karyoplot, transcript = transcript, 
                       coding.exons = coding.exons, non.coding.exons = non.coding.exons, 
                       r0=transcript.y0, r1=transcript.y1, col=tcol)
      }
    }
      
      #t.exons <- unlist(all.exons[transcript$tx_id])
      
      
      transcript.r0 <- bin.height*(bin-1)
      transcript.r1 <- transcript.r0 + transcript.height
      
      #kpRect(kp, transcript, y0=0, y1=1, r0=transcript.r0, r1=transcript.r1, col=col)  
      
      plotTranscript(kp, transcript = transcript, coding.exons = coding.exons, non.coding.exons = non.coding.exons, 
                     strand = as.character(strand(transcript))[1],
                     transcript.name.col="#666666", non.coding.exons.height = 0.5, strand.marks = FALSE, 
                     mark.height = NULL, mark.width = NULL, mark.distance = 4, transcript.name = transcript$tx_name, transcript.name.position = "left", transcript.name.cex = 1, transcript.name.offset = 0, coding.exons.col = tcol, non.coding.exons.col = tcol, marks.col = tcol, introns.col=tcol,
                     r0=transcript.r0, r1=transcript.r1, data.panel=1)
      
    }
    
    
    
    
    karyoplot <- plotKaryotype(genome="hg19", zoom=zoom)
    kpAddBaseNumbers(karyoplot, tick.dist = 100e3)
    
    
    
    
    
    
    if(plot.transcripts.structure==FALSE) { #Plot the transcripts as boxes
      
    }
  }
  
  
 
  
  
  karyoplot <- plotKaryotype(genome="hg19", zoom=zoom)
  #kpAddBaseNumbers(kp, tick.dist = 5e3)
  
  
  
  
  
  
  invisible(karyoplot)
}



#TODO: Make it exportable so users can build the structures from it and modify them as needed. (Changing the gene names?)
makeDataFromTxDb <- function(karyoplot, txdb, plot.transcripts, plot.transcripts.structure) {
  res <- list()
  
  #get the genes
  all.genes <- genes(txdb)
  res[["genes"]] <- subsetByOverlaps(all.genes, karyoplot$plot.region)
  
  if(plot.transcripts==TRUE) {
    #Extract the transcripts of each gene
    res[["transcripts"]] <- as.list(transcriptsBy(txdb, by="gene")[names(res$genes)])
  
    if(plot.transcripts.structure==TRUE) {
      #extract the exons of each transcript
      res[["coding.exons"]] <- list()
      res[["non.coding.exons"]] <- list()
      #for every transcript, get the coding exons and compute the non coding ones
      all.trans <- as.character(Reduce(c, res$transcripts)$tx_id)
      for(tid in all.trans) {
        res$coding.exons[[tid]] <- cds(txdb, filter=list(tx_id=tid))
        complete.exons <- exons(txdb, filter=list(tx_id=tid))
        res$non.coding.exons[[tid]] <- subtractRegions(complete.exons, res$coding.exons[[tid]])
      }
    }
  }
  return(res)
}

#TODO
isValidData <- function(...) {
  return(TRUE)
}


#plot names of genes or transcripts (given the rectangle they are labelling)
plotNames <- function(karyoplot, data, y0, y1, labels, position, col, cex, r0, r1, data.panel) {
  if(position=="left") {
    message("y=", (y1-y0)/2)
    kpText(karyoplot, chr=as.character(seqnames(data)), x=start(data), y=y0+(y1-y0)/2,  labels=labels, pos=2, cex=cex, col=col, r0=r0, r=r1, data.panel)
  }
  if(position=="right") {
    kpText(karyoplot, chr=as.character(seqnames(data)), x=end(data), y=y0+(y1-y0)/2,  labels=labels, pos=4, cex=cex, col=col, r0=r0, r=r1, data.panel)
  }
  if(position=="top") {
    kpText(karyoplot, chr=as.character(seqnames(data)), x=start(data)+width(data)/2, y=y1,  labels=labels, pos=3, cex=cex, col=col, r0=r0, r=r1, data.panel)
  }
  if(position=="bottom") {
    kpText(karyoplot, chr=as.character(seqnames(data)), x=start(data)+width(data)/2, y=y0,  labels=labels, pos=1, cex=cex, col=col, r0=r0, r=r1, data.panel)
  }
}




