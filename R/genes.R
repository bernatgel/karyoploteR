################################################################################
#  The functions in this file are utility functions to build the GenesData
#  structure needed by kpPlotGenes from different sources and to manage
#  the conversion of gene identifiers to gene symbols
################################################################################

#Dataframe tying the available OrgDb objects to the different identifiers of 
#the organism they annotate
#Data extracted from https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
.OrganismToOrgDb <- data.frame(
  taxonomyId=c(9606, 10090),
  organism=c("Homo sapiens", "Mus musculus"),
  genome=c("hg", "mm"),
  package=c("org.Hs.eg.db", "org.Mm.eg.db"),
  stringsAsFactors=FALSE
)



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
#' @importFrom GenomicFeatures genes exons transcriptsBy cds organism
#' @importFrom IRanges subsetByOverlaps
#' @importFrom AnnotationDbi taxonomyId select
#' @importFrom GenomeInfoDb genome
#'

makeGenesDataFromTxDb <- function(karyoplot, txdb, plot.transcripts=TRUE, plot.transcripts.structure=TRUE) {
  res <- list()
  
  #get the metadata
  res$metadata <- list()
  res$metadata$organism <- GenomicFeatures::organism(txdb)
  res$metadata$taxonomyId <- AnnotationDbi::taxonomyId(txdb)
  res$metadata$genome <- setNames(GenomeInfoDb::genome(txdb)[1], NULL)
  
  
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
  class(res) <- "GenesData"
  return(res)
}



#' addGeneNames
#' 
#' @description 
#' Adds the gene names (defaults to symbols) to a GenesData object to 
#' be used by kpPlotGenes
#' 
#' 
#' @details 
#' This function takes a valid data object and uses an OrgDb object to 
#' find the gene names (symbols by default) and add them. Names are added 
#' as a column named \code{names} to the \code{genes} element of GenesData 
#' and they replace anything that was present there before.
#' If no \code{ObjDb} object is given, the function will try to identify
#' the correct organism using the data in \code{GenesData$metadata} and
#' select the OrgDb object if available. If it cannot identify the organism 
#' or there's no valid OrgDb for that organism it will fail with an error.
#' Internally, the function uses a call to \code{AnnotationDbi::select} on
#' the OrgDb. It is possible to specify the keys and keytypes as well as
#' the column we want to use as names (defaults to SYMBOL for gene symbols).
#' 
#' @usage addGeneNames(genes.data, orgDb="auto", keys=NULL, keytype="ENTREZID", names="SYMBOL")
#' 
#' @param genes.data (GenesData object) A valid genes.dat object like the ones obtained by \code{\link{makeGenesDataFromTxDb}}
#' @param orgDb The orgDb object to use to extract the gene symbols. If "auto" the function will try to determine automatically the correct organism. See available obects at https://bioconductor.org/packages/release/BiocViews.html#___OrgDb (defaults to "auto")
#' @param keys (character vector ) The keys to be used in the internal select statement to get the names. If NULL, the first column of \code{mcols(GenesData$genes)} will be used. (defaults to NULL)
#' @param keytype (character) The keytype used in the internal select statement. (defaults to "ENTREZID", that is, gene_id)
#' @param names The column to extract from orgDb to use as gene names. (deafults to "SYMBOL")
#' 
#' 
#' @return
#' The original GenesData object with one additional column named "names" in 
#' \code{GenesData$genes$names}.
#' 
#'  
#' @seealso \code{\link{kpPlotGenes}}, \code{\link{makeGenesDataFromTxDb}}
#' 
#' @examples
#'  
#' 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' zoom <- toGRanges("chr17:29e6-30e6")
#' kp <- plotKaryotype(genome="hg19", zoom=zoom)
#' genes.data <- makeGenesDataFromTxDb(kp, TxDb.Hsapiens.UCSC.hg19.knownGene, FALSE, FALSE)
#' genes.data <- addGeneNames(genes.data)
#' kpPlotGenes(kp, data=genes.data, r1=0.5, plot.transcripts=FALSE, gene.name.position = "left")
#'  
#'   
#' @export addGeneNames
#' 
addGeneNames <- function(genes.data, orgDb="auto", keys=NULL, keytype="ENTREZID", names="SYMBOL") {
  if(!methods::is(genes.data, "GenesData")) stop("genes.data must be a valid object of the GenesData class")
  #TODO: check it's a valid genes.data
  
  if(is.character(orgDb) & orgDb=="auto") {
    org <- NULL
    if(!is.null(genes.data$metadata$taxonomyId)) {
      org <- .OrganismToOrgDb[.OrganismToOrgDb$taxonomyId==genes.data$metadata$taxonomyId,]
    }
    if(is.null(org) && !is.null(genes.data$metadata$genome)) {
      
      #Use the genome name (prefix) to get the correct org line
    }
    
    
    if(!is.null(org)) {
      loaded <- require(org$package, character.only = TRUE)
      if(loaded==TRUE) orgDb <- get(org$package)
    } else {
      stop("It was not possible to identify the organism in order to select the appropiate OrgDb object. Please provide one manually")
    }
  }
  if(!methods::is(orgDb, "OrgDb")) stop("orgDb must be either a valid OrgDb object or 'auto'")
 
  #If we are here, we have a valid orgDb 
  if(is.null(keys)) keys <- mcols(genes.data$genes)[,1]
  
  ss <- suppressMessages(AnnotationDbi::select(orgDb, keys=keys, keytype=keytype, columns=names))
  
  genes.data$genes$name <- ss[,2]

  return(genes.data)  
}



#' mergeTranscripts
#' 
#' @description 
#' Merges the transcripts of each gene and creates one transcript per gene with
#' all exons and UTR regions combined
#' 
#' 
#' @details 
#' This function takes a valid data object and merges all transcripts from 
#' each gene into a single transcript. This is useful to reduce the plot
#' complexity while keeping partial information on transcript structure.#' 
#' In this transcript, any region that is a coding exon in any transcript,
#' will be an exon, any region that is a non-coding exon in any transcript
#' and is not an exon in any transcript, will be a non-coding exon. 
#' Anything between coding and non-coding exons will be introns.
#' 
#' @usage mergeTranscripts(genes.data)
#' 
#' @param genes.data (GenesData object) A valid genes.dat object like the ones obtained by \code{\link{makeGenesDataFromTxDb}}
#' 
#' @return
#' The original GenesData object with a single transcript per gene 
#' \code{GenesData$genes$names}.
#' 
#' @seealso \code{\link{kpPlotGenes}}, \code{\link{makeGenesDataFromTxDb}}
#' 
#' @examples
#'  
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 
#' zoom <- toGRanges("chr17:29.4e6-29.8e6")
#' kp <- plotKaryotype(genome="hg19", zoom=zoom)
#' genes.data <- makeGenesDataFromTxDb(kp, TxDb.Hsapiens.UCSC.hg19.knownGene)
#' genes.data <- addGeneNames(genes.data)
#' kpPlotGenes(kp, data=genes.data, r1=0.5, plot.transcripts=TRUE, gene.name.position = "left")
#' genes.data.merged <- mergeTranscripts(genes.data)
#' kpPlotGenes(kp, data=genes.data.merged, r0=0.6, r1=0.8, plot.transcripts=TRUE, gene.name.position = "left")
#'  
#' @export mergeTranscripts
mergeTranscripts <- function(genes.data) { 
  if(!methods::is(genes.data, "GenesData")) stop("genes.data must be a valid object of the GenesData class")
  #TODO: check it's a valid genes.data
  
  if(is.null(genes.data$transcripts)) stop("genes.data must have transcript information. If created with makeGenesDataFromTxDb,  plot.transcripts must be set to TRUE")
  if(is.null(genes.data$coding.exons)) stop("genes.data must have exons information. If created with makeGenesDataFromTxDb,  plot.transcripts.structure must be set to TRUE")

  merged.data <- list()
  class(merged.data) <- "GenesData"
  merged.data$genes <- genes.data$genes
  merged.data$transcripts <- list()
  merged.data$coding.exons <- list()
  merged.data$non.coding.exons <- list()
  
  for(ngene in seq_len(length(genes.data$genes))) {
    g <- genes.data$genes[ngene]
  
    #Get the transcripts and get the total region overlapping any of them
    tx <- genes.data$transcripts[[names(g)]]
    
    merged.tx <- regioneR::joinRegions(tx, 0)
    merged.tx_id <- paste0(names(g), ".merged")
    merged.tx$tx_id <- merged.tx_id
    
    merged.data$transcripts[[names(g)]] <- merged.tx
    
    #Get all coding exons and merge them
    cod.ex <- do.call("c", lapply(tx$tx_id, function(tx_id) return(genes.data$coding.exons[[as.character(tx_id)]])))
    cod.ex <- regioneR::joinRegions(cod.ex, 0)
    
    merged.data$coding.exons[[merged.tx_id]] <- cod.ex
    
    #Get all non-coding exons and merge them
    ncod.ex <- do.call("c", lapply(tx$tx_id, function(tx_id) return(genes.data$non.coding.exons[[as.character(tx_id)]])))
    ncod.ex <- regioneR::joinRegions(ncod.ex, 0)
    
    merged.data$non.coding.exons[[merged.tx_id]] <- ncod.ex
  }
  return(merged.data) 
}



#TODO
isValidData <- function(...) {
  return(TRUE)
}





getGeneNames <- function(genes, gene.names) {
  if(!is.null(gene.names)) {
    return(naToEmptyChar(as.character(gene.names[names(genes)])))
  } else {
    if("name" %in% names(mcols(genes))) {
      return(naToEmptyChar(genes$name))
    } else {
      if(length(mcols(genes))>0) {
        return(naToEmptyChar(as.character(mcols(genes)[,1])))
      } else {
        return(naToEmptyChar(as.character(names(genes))))
      }
    }
  }
}

getTranscriptNames <- function(transcripts, transcript.names) {
  if(!is.null(transcript.names)) {
    return(as.character(transcript.names[names(transcripts)]))
  } else {
    return(as.character(names(transcripts)))
  }
}
