################################################################################
#  The functions in this file are utility functions to build the GenesData
#  structure needed by kpPlotGenes from different sources and to manage
#  the conversion of gene identifiers to gene symbols
################################################################################

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
#' @usage makeGenesDataFromTxDb(karyoplot, txdb, plot.transcripts=TRUE, plot.transcripts.structure=TRUE)
#' 
#' @param karyoplot  (karyoplot object) A valid karyoplot object created by a call to \code{\link{plotKaryotype}}
#' @param txdb (a TxDb object) The transcript database object from which the data will be extracted.
#' @param plot.transcripts (boolean) TRUE if transcripts are needed in addition to the genes. (defaults to TRUE)
#' @param plot.transcripts.structure (boolean) TRUE if the coding and non-coding exons are needed in addition to the genes and transcripts. (defaults to TRUE)
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

makeGenesDataFromTxDb <- function(karyoplot, txdb, plot.transcripts=TRUE, plot.transcripts.structure=TRUE) {
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
