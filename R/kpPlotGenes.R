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


kpPlotGenes <- function(karyoplot, data, gene.margin=0.3, gene.col=NULL, gene.border.col=NULL,
                        add.gene.names=TRUE, gene.names=NULL, gene.name.position="top", gene.name.cex=1, gene.name.col=NULL,
                        plot.transcripts=TRUE, transcript.margin=0.5, transcript.col=NULL, transcript.border.col=NULL,
                        add.transcript.names=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=0.6, transcript.name.col=NULL,
                        plot.transcripts.structure=TRUE,
                        non.coding.exons.height=0.5, 
                        add.strand.marks=TRUE, mark.height=0.20, mark.width=1, mark.distance=4,
                        coding.exons.col=NULL, coding.exons.border.col=NULL, 
                        non.coding.exons.col=NULL, non.coding.exons.border.col=NULL, 
                        introns.col=NULL, marks.col=NULL,
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
    data <- tryCatch(makeGenesDataFromTxDb(karyoplot, data, plot.transcripts, plot.transcripts.structure),
                     error=function(e) {stop("Error: There was an error extracting the information from the TxDb object. ", e)})
  }
  
  #if there's nothing to plot, return
  if(length(data$genes)==0) {
    #add an empty "latest.plot" information
    #TODO
    #and return
    invisible(karyoplot)
  }
  
  if(length(data$genes)>20 & plot.transcripts==TRUE & plot.transcripts.structure==TRUE) {
    message("NOTE: Plotting many genes with detailed transcript structure may take a long time. You can set 'plot.transcripts' and 'plot.transcripts.structure' to FALSE to speed up the process reducing the detail level.")
  }
  
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  total.height <- 1
  
  #Set the colors for the unspecified entities
  if(is.null(border)) {border <- col}
  if(is.null(gene.col)) { gene.col<- col}
  if(is.null(gene.border.col)) { gene.border.col<- border}
  if(is.null(gene.name.col)) { gene.name.col<- col}
  if(is.null(transcript.col)) { transcript.col<- col}
  if(is.null(transcript.border.col)) { transcript.border.col<- border}
  if(is.null(transcript.name.col)) { transcript.name.col<- col}
  if(is.null(coding.exons.col)) { coding.exons.col<- col}
  if(is.null(coding.exons.border.col)) { coding.exons.border.col<- border} 
  if(is.null(non.coding.exons.col)) { non.coding.exons.col<- col} 
  if(is.null(non.coding.exons.border.col)) { non.coding.exons.border.col<- border} 
  if(is.null(introns.col)) { introns.col<- col}
  if(is.null(marks.col)) { marks.col <- col}

  
  #All parameters processed. Start working
  gene.pos <- list()
  transcript.pos <- list()
  
  
  if(plot.transcripts==FALSE) {
    #Plot only the genes, with one rectangle per gene. Automatically position genes so they do not overlap
    if(avoid.overlapping==TRUE) {
      genes.for.coverage <- data$genes
      strand(genes.for.coverage) <- "*"
      bins <- disjointBins(genes.for.coverage)
      num.layers <- max(bins)
      layer.height <- total.height/num.layers
      gene.height <- layer.height/(1+gene.margin)
      y0 <- layer.height * (bins-1)
      y1 <- y0 + gene.height
    } else {
      y0 <- 0
      y1 <- 1 - gene.margin
    }
    kpRect(karyoplot, data=data$genes, y0=y0, y1=y1, col=gene.col, border=gene.border.col, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping, ...)
    if(add.gene.names==TRUE) {
      gene.labels <- getGeneNames(genes=data$genes, gene.names=gene.names)
      kpPlotNames(karyoplot, data=data$genes, y0=y0, y1=y1, labels=gene.labels, position=gene.name.position, col=gene.name.col, cex=gene.name.cex, r0=r0, r1=r1, clipping=clipping, data.panel=data.panel, ...)
    }
    
  } else { #if plot.transcripts==TRUE
    #Compute the position of each transcript (and the height dedicated to each 
    # gene depending on the numbre of transcripts it has) and then, plot them
    
    if(avoid.overlapping==FALSE) {
      #All transcripts span the whole vertical space
      all.transcript.names <- as.character(unlist(GRangesList(data$transcripts))$tx_id)
      for(t in all.transcript.names) {
        transcript.pos[[t]] <- list(y0=0, y1=1)
      }
    } else { #Get the positions of the transcripts so they don't overlap
      #TODO: Add the transcript name length to avoid overlapping?
      
      #Treat all genes as if they have the same length as the gene, so each genes has its own "rectangular space" preserved
      num.transcripts.per.gene <- Map(length, as.list(data$transcripts))
      genes.for.coverage <- Reduce(c, Map(rep, as.list(data$genes), num.transcripts.per.gene))
      strand(genes.for.coverage) <- "*"
      max.coverage <- max(max(coverage(genes.for.coverage)))
      
      bin.height <- total.height/max.coverage
      transcript.height <- bin.height/(1+transcript.margin)
      
      #Position them using the disjoint binning functionality of Bioconductor
      bins <- disjointBins(genes.for.coverage)
      names(bins) <- names(genes.for.coverage)
      
      
      #TODO: This could be vectorized with a tapply per gene in the bins object
      # and inside, an apply to each transcript to get its position
      #Initialize the transcript counter
      transcript.in.gene <- 1
      last.gene <- ""
      for(nt in seq_len(length(bins))) {
        bin <- bins[nt]
        #if we are starting a new gene restart the transcript counter, 
        #else, move to the next transcript in the gene
        if(last.gene != names(bin)) { 
          transcript.in.gene <- 1
          last.gene <- names(bin)
        } else {
          transcript.in.gene <- transcript.in.gene + 1
        }
        
        #And compute the y0 and y1 of that transcript
        transcript <- data$transcripts[[names(bin)]][transcript.in.gene]
        names(transcript) <- transcript$tx_id #Does it work with all txdbs??
        
        transcript.y0 <- bin.height*(bin-1)
        transcript.y1 <- transcript.y0 + transcript.height
        
        transcript.pos[[as.character(names(transcript))]] <- list(y0=setNames(transcript.y0, NULL), y1=setNames(transcript.y1, NULL))
        
      }
      
    }
    
    #Once we've got the position of transcripts, compute the position of the genes  
    for(g in as.character(names(data$genes))) {
      #get the vertical position of all transcripts of that gene
      tpos <- data.frame(do.call(rbind, transcript.pos[as.character(data$transcripts[[g]]$tx_id)]))
      gene.pos[[g]] <- list(y0=min(unlist(tpos[,1])), y1=max(unlist(tpos[,2])))
    }
    
    #All positioning data ready. Plot the transcripts
    
    #Now, plot the transcripts in their positions
    detail.level <- ifelse(plot.transcripts.structure==TRUE, 2, 1)
      
    #Join all transcripts in a single GRanges as expected by kpPlotTranscripts
    transcripts <- unlist(GRangesList(data$transcripts))
    names(transcripts) <- transcripts$tx_id
    
    #Compute the y0 and y1 of all transcripts
    y0 <- unlist(do.call(rbind, transcript.pos[names(transcripts)])[,1])
    y1 <- unlist(do.call(rbind, transcript.pos[names(transcripts)])[,2])
    
    kpPlotTranscripts(karyoplot, 
                      data=list(transcripts=transcripts, coding.exons=data$coding.exons, non.coding.exons=data$non.coding.exons), 
                      add.transcript.names=add.transcript.names, transcript.names = transcript.names,
                      y0 = y0, y1=y1, detail.level = detail.level, 
                      r0=r0, r1=r1, 
                      non.coding.exons.height = non.coding.exons.height, 
                      add.strand.marks = add.strand.marks, mark.height = mark.height, mark.width = mark.width,
                      mark.distance = mark.distance, marks.col = marks.col, 
                      transcript.name.position = transcript.name.position, 
                      transcript.name.cex = transcript.name.cex, 
                      transcript.name.col = transcript.name.col, 
                      coding.exons.col = coding.exons.col, 
                      coding.exons.border.col = coding.exons.border.col,
                      non.coding.exons.col = non.coding.exons.col, 
                      non.coding.exons.border.col = non.coding.exons.border.col, 
                      introns.col = introns.col,
                      col=transcript.col,
                      border=transcript.border.col,
                      data.panel = data.panel, clipping = clipping,
                      ymin=0, ymax=1, ...) 
                        
    #Finally, plot the gene names if needed
    if(add.gene.names) {
      gene.labels <- getGeneNames(data$genes, gene.names)
      y0 <- unlist(do.call(rbind, gene.pos[as.character(names(data$genes))])[,1])
      y1 <- unlist(do.call(rbind, gene.pos[as.character(names(data$genes))])[,2])
      kpPlotNames(karyoplot, data=data$genes, y0=y0, y1=y1, labels=gene.labels, col=gene.name.col, cex=gene.name.cex, position = gene.name.position, r0 = r0, r1=r1, clipping=clipping, data.panel=data.panel)
    }
  } #END if plot.transcripts==TRUE
  
  #Store the gene and transcript position in the karyoplot object
  karyoplot$latest.plot <- list(funct="kpPlotGenes", 
                                computed.values=list(gene.vertical.position=gene.pos, 
                                                     transcript.vertical.position=transcript.pos)
  )

  invisible(karyoplot)
}


#TODO: Make it exportable so users can build the structures from it and modify them as needed. (Changing the gene names?)
makeGenesDataFromTxDb <- function(karyoplot, txdb, plot.transcripts, plot.transcripts.structure) {
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





getGeneNames <- function(genes, gene.names) {
  if(!is.null(gene.names)) {
    return(as.character(gene.names[names(genes)]))
  } else {
    return(as.character(names(genes)))
  }
}

getTranscriptNames <- function(transcripts, transcript.names) {
  if(!is.null(transcript.names)) {
    return(as.character(transcript.names[names(transcripts)]))
  } else {
    return(as.character(names(transcript)))
  }
}
