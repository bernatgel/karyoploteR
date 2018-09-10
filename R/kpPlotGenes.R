#' kpPlotGenes
#' 
#' @description 
#' 
#' Plot genes and transcripts in the genome. Can get the genes and trancripts information from TxDb or from custom objects.
#' 
#' @details 
#'
#'  This is one of karyoploteR's higher level functions. It takes a transcript 
#'  database (\code{TxDb}) object or a custom object with a specific structure 
#'  and plots the genes along the genome. It's possible to plot genes as 
#'  a whole using rectangle for each gene or to plot each transcript 
#'  independently. If transcripts are drawn, it's possible to plot them as
#'  single boxes or to plot the detailed structure, differentiating 
#'  coding exons, non-coding exons and introns. Transcripts may have, in 
#'  addition, little arrows to mark the trascript strand (plus or minus). These
#'  strand marks are plotted in introns if the transcript structure is shown 
#'  or on the whole transcript length if transcripts are plotted as boxes. 
#'  Finally, it's possible to add labels to genes and transcripts. By default 
#'  genes and transcripts identifiers in the input data structure will be used
#'  as labels, but it's possible to provide named character vectors to be 
#'  used as dictionaries to change id's to better names.
#'    
#'  The genes and transcripts representations are customizable. It's possible
#'  to change the colors of the different elements individually (i.e. to 
#'  have red coding exons, blue non-coding exons and green introns); it's 
#'  possible to change the relative height of the non-coding exons and to 
#'  change the slate and density of the strand marks.
#'  
#'  The data stating the positions of genes, transcripts and exons in the 
#'  genome and their relations (which transcripts belong to which genes) 
#'  can be given as a standard transcript database (\code{TxDb}) object or
#'  as a custom list with the following elements: \code{genes},
#'  \code{transcripts}, \code{coding.exons} and \code{non.coding.exons}.
#'  
#'  
#' @note Plotting transcripts, specially plotting their structure might get 
#' quite slow in comparison to the usual speed of plotting in karyoploteR. It
#' is not advised to plot genes and transcripts on the whole genome or in 
#' large regions of it. These functions have been designed to work with
#' zoomed in karyoplots.
#'
#' @usage kpPlotGenes(karyoplot, data, gene.margin=0.3, gene.col=NULL, gene.border.col=NULL,
#'                        add.gene.names=TRUE, gene.names=NULL, gene.name.position="top", gene.name.cex=1, gene.name.col=NULL,
#'                        plot.transcripts=TRUE, transcript.margin=0.5, transcript.col=NULL, transcript.border.col=NULL,
#'                        add.transcript.names=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=0.6, transcript.name.col=NULL,
#'                        plot.transcripts.structure=TRUE,
#'                        non.coding.exons.height=0.5, 
#'                        add.strand.marks=TRUE, mark.height=0.20, mark.width=1, mark.distance=4,
#'                        coding.exons.col=NULL, coding.exons.border.col=NULL, 
#'                        non.coding.exons.col=NULL, non.coding.exons.border.col=NULL, 
#'                        introns.col=NULL, marks.col=NULL,
#'                        data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col="black", 
#'                        border=NULL, avoid.overlapping=TRUE, clipping=TRUE, ...)
#' 
#' 
#' @param karyoplot (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data  (a \code{TxDb} object or a list with the required elements) A \code{TxDb} object with information on genes, transcripts and exons and their position on the genome.
#' @param gene.margin  (numeric) If whole genes (as oposed to transcripts) are plotted (\code{plot.transcripts=FALSE}), the vertical margin between overlapping genes. The value is with respect to the gene height. A value of 0.5 will create a space above the genes with half the height of the genes themselves. (defaults to 0.3)
#' @param gene.col  (color) If whole genes (as oposed to transcripts) are plotted (\code{plot.transcripts=FALSE}), the color used to fill the rectangles representing the genes. If NULL, the value of \code{col} will be used. (Defaults to NULL)
#' @param gene.border.col  (color) If whole genes (as oposed to transcripts) are plotted (\code{plot.transcripts=FALSE}), the color used in the border of the rectangles representing the genes. If NULL, the value of \code{col} will be used. If NA, no border is drawn. (Defaults to NULL)
#' @param add.gene.names  (boolean) Whether to add the names of the genes to the plot.
#' @param gene.names  (named character vector) A named character vector with the labels of the genes. If not NULL, it will be used as a dictionary, so gene ids should be names and desired labels the values. If NULL, the gene ids will be used as labels. (defaults to null) 
#' @param gene.name.position (character) The position of the gene name text relative to the rectangle. Can be "left", "right", "top", "bottom" or "center". (Defaults to "top")
#' @param gene.name.cex (numeric) The cex value to plot the gene names (defaults to 1)
#' @param gene.name.col (color) The color of the gene labels. If NULL, it will use col. (defaults to NULL)
#' @param plot.transcripts (boolean) Whether to plot the individual transcripts (TRUE) or a single rectangle to represent the whole gene (FALSE). (Defaults to TRUE)
#' @param transcript.margin (numeric) If transcripts are plotted (\code{plot.transcripts=TRUE}), the vertical margin between overlapping transcripts. The value is with respect to the transcript height. A value of 0.5 will create a space above the transcripts with half the height of the transcripts themselves. (defaults to 0.5)
#' @param transcript.col (color) If transcripts are plotted (\code{plot.transcripts=TRUE}), the color used to fill the rectangles representing the transcripts. If NULL, the value of \code{col} will be used. (Defaults to NULL)  
#' @param transcript.border.col (color) If transcripts are plotted (\code{plot.transcripts=TRUE}), the color used in the border the rectangles representing the transcripts. If NULL, the value of \code{col} will be used. If NA, no border will be drawn. (Defaults to NULL)  
#' @param add.transcript.names (boolean) Whether to add the names of the tramscripts to the plot.
#' @param transcript.names (named character) A named character vector with the labels of the transcripts. If not null, it will be used as a dictionary, so transcript ids should be names and desired labels the values. If NULL, the transcript ids will be used as labels. (defaults to null)
#' @param transcript.name.position (character) The position of the transcript name text relative to the rectangle. Can be "left", "right", "top", "bottom" or "center". (Defaults to "left")
#' @param transcript.name.cex (numeric) The cex value to plot the transcript names (defaults to 0.6)
#' @param transcript.name.col (color) The color of the transcript labels. If NULL, it will use col. (defaults to NULL)
#' @param plot.transcripts.structure (boolean) Whether to draw the transcripts as single rectangles (FALSE) or to show the complete transcript structure (introns and exons) (TRUE). (Defaults to TRUE)
#' @param non.coding.exons.height  (numeric) The height of the non.coding exons relative to the transcript height. For example, if 0.5, non-coding exons will have a height half the size of the coding ones. (default 0.5) 
#' @param add.strand.marks  (boolean) Whether strand marks should be plotted or not. Strand marks are small arrows along the introns (or whole transcripts if plot.transcript.structure=FALSE). (defaults to TRUE)
#' @param mark.height (numeric) The height of the strand marks in "coding exons heights", that is, if mark.height is 0.5, the mark will have a height of half the height of an exon. (defaults to 0.2)
#' @param mark.width (numeric) The width of the strand marks, in mark heights. mark.width=1 will produce arrow heads with a slope pf 45 degrees. A value higher than 1 will produce smaller angles and a value below 1 larger angles with more vertical lines. (defaults to 1, 45 degrees)
#' @param mark.distance   (numeric) The distance between marks, in mark widths. A distance of 2, will add a space of 2*mark.width between consecutive marks. (defaults to 4)
#' @param coding.exons.col  (color) The fill color of the rectangles representing the coding exons. If NULL, it will use col. (defaults to NULL)
#' @param coding.exons.border.col  (color) The color of the border of the coding exons. If NULL, it will use border. (defaults to NULL)
#' @param non.coding.exons.col  (color) The fill color of the rectangles representing the non-coding exons. If NULL, it will use col. (defaults to NULL)
#' @param non.coding.exons.border.col  (color) The color of the border of the non-coding exons. If NULL, it will use border. (defaults to NULL)
#' @param introns.col (color) The color of the lines representing the introns. If NULL, it will use col. (defaults to NULL) 
#' @param marks.col   (color) The color of the arrows representing the strand. If NULL, it will use col. (defaults to NULL) 
#' @param data.panel (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0  (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL) 
#' @param r1  (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL) 
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track min(autotrack[[1]]) to track max(autotrack[[1]]). The margin is specified as the part of a track, by default 0.05, 5% of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param col  (color) The color of the genes, transcripts and labels. It is possible to specify different colors for each element class (transcript names, exons, strand marks...). All elements with no explicit color will be plotted using col. (Defaults to "black") 
#' @param border  (color) The color of the border of rectangles representing genes, transcripts and exons. Every element class may have its own specific color using the appropiate parameters. The ones with no explicit color will use border. At the same time, if border is NULL, it will default to col. (Defaults to NULL)
#' @param avoid.overlapping (boolean) If two or more regions (genes, transcripts) overlap in the genome (even partially), they can be drawn in two different layers (TRUE) or in the same layer, with superposing representations (FALSE). (Defaults to TRUE, draw non-overlapping elements)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...  The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
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
#'  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'  
#'  zoom <- toGRanges("chr2", 47986268, 48147403)
#'  gene.names <- c("2956"="MSH6", "80204"="FBXO11")
#'  
#'  kp <- plotKaryotype(genome="hg19", zoom=zoom)
#'  kpPlotGenes(kp, data=txdb, add.transcript.names = FALSE, gene.names=gene.names, r1=0.6)
#' 
#'  kp <- plotKaryotype(genome="hg19", zoom=zoom)
#'  kpPlotGenes(kp, data=txdb, plot.transcripts=FALSE, gene.names=gene.names, r1=0.6)
#'
#'  kp <- plotKaryotype(genome="hg19", zoom=zoom)
#'  kpPlotGenes(kp, data=txdb, plot.transcripts.structure=FALSE, add.transcript.names=FALSE, gene.names=gene.names, r1=0.8, col="blue", marks.col="white", gene.name.col="black")
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
                        data.panel=1, r0=NULL, r1=NULL, autotrack=NULL, col="black", 
                        border=NULL, avoid.overlapping=TRUE, clipping=TRUE, ...) {
  
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
    kpRect(karyoplot, data=data$genes, y0=y0, y1=y1, col=gene.col, border=gene.border.col, r0=r0, r1=r1, autotrack=autotrack, data.panel=data.panel, clipping=clipping, ...)
    if(add.gene.names==TRUE) {
      gene.labels <- getGeneNames(genes=data$genes, gene.names=gene.names)
      kpPlotNames(karyoplot, data=data$genes, y0=y0, y1=y1, labels=gene.labels, position=gene.name.position, col=gene.name.col, cex=gene.name.cex, r0=r0, r1=r1, autotrack=autotrack, clipping=clipping, data.panel=data.panel, ...)
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
      num.transcripts.per.gene <- lengths(data$transcripts)
      genes.for.coverage <- rep(data$genes, num.transcripts.per.gene)
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
                      r0=r0, r1=r1, autotrack=autotrack,
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
      kpPlotNames(karyoplot, data=data$genes, y0=y0, y1=y1, labels=gene.labels, col=gene.name.col, cex=gene.name.cex, position = gene.name.position, r0 = r0, r1=r1, autotrack=autotrack, clipping=clipping, data.panel=data.panel)
    }
  } #END if plot.transcripts==TRUE
  
  #Store the gene and transcript position in the karyoplot object
  karyoplot$latest.plot <- list(funct="kpPlotGenes", 
                                computed.values=list(gene.vertical.position=gene.pos, 
                                                     transcript.vertical.position=transcript.pos)
  )

  invisible(karyoplot)
}


#' makeGenesDataFromTxDb
#' 
#' @description 
#' 
#' This is a utility function that transforms a TxDb object into a custom 
#' object valid as input for \code{\link{kpPlotGenes}}.
#' 
#' @details 
#'  
#' This function creates a valid data object for \code{\link{kpPlotGenes}} 
#' starting from a TxDb object. The resulting object will contain only the
#' genes and transcripts ovelapping the plot region of the given Karyoplot 
#' object. 
#' 
#' @note \code{\link{kpPlotGenes}} accepts TxDb objects directly. This 
#' function is only expected to be used when the user want to manipulate the 
#' results somehow (i.e. removing some of the genes).
#' 
#' @usage makeGenesDataFromTxDb(karyoplot, txdb, plot.transcripts, plot.transcripts.structure)
#' 
#' @param karyoplot  (karyoplot object) A valid karyoplot object created by a call to \code{\link{plotKaryotype}}
#' @param txdb (a TxDb object) The transcript database object from which the data will be extracted.
#' @param plot.transcripts (boolean) TRUE if transcripts are needed in addition to the genes.
#' @param plot.transcripts.structure (boolean) TRUE if the coding and non-coding exons are needed in addition to the genes and transcripts.
#' 
#' @return
#' 
#' Returns a list with at least one element called \code{genes}, a 
#' \code{GRanges} with all genes overlapping karyoplot. If 
#' \code{plot.transcripts} is TRUE, the returned list will have 
#' a \code{transcript} element, a list of \code{GRanges} objects, one per gene
#' (named with the gene ids),
#' with the transcripts of that gene. If \code{plot.transcripts.structure} is
#' TRUE, two more elements are present: \code{coding.exons} and 
#' \code{non.coding.exons}, each a list with one element per trascript 
#' (named with the transcript id), and each element the coding or non-coding
#' exons of that transcript.
#' 
#'  
#' @seealso \code{\link{kpPlotGenes}}
#' 
#' @examples
#'  
#' 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' zoom <- toGRanges("chr17", 32.6e6, 33.2e6)
#' kp <- plotKaryotype(genome="hg19", zoom=zoom)
#' 
#' genes.data <- makeGenesDataFromTxDb(kp, TxDb.Hsapiens.UCSC.hg19.knownGene, TRUE, TRUE)
#'   
#' @export makeGenesDataFromTxDb
#' 
#' @importFrom GenomicFeatures genes exons transcriptsBy cds
#' @importFrom IRanges subsetByOverlaps
#' 


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
    return(as.character(names(transcripts)))
  }
}
