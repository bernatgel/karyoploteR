#' getCytobands
#' 
#' @description 
#' 
#' Get the cytobands of the specified genome.
#' 
#' @details 
#'  
#' It returns \code{GRanges} object with the cytobands of the specified genome. 
#' The cytobands for some organisms and genome versions have been pre-downloaded from UCSC
#' and included in the \code{karyoploteR} package. For any other genome, \code{getCytobands}
#' will use \code{rtracklayer} to try to fetch the \code{cytoBandIdeo} table from UCSC. If for
#' some reason it is not possible to retrieve the cytobands, it will return an empty \code{GRanges}
#' object. Setting the parameter \code{use.cache} to \code{FALSE}, the data included in the 
#' package will be ignored and the cytobands will be downloaded from UCSC.
#' 
#' The genomes (and versions) with pre-downloaded cytobands are: hg19, hg38, mm9, mm10, rn5, rn6,
#' danRer10, dm6, ce6 and sacCer3.
#'    
#' 
#' @usage getCytobands(genome="hg19", use.cache=TRUE)
#' 
#' @param genome   (character or other) specifies a genome using the UCSC genome name. Defaults to "hg19". If it's not a \code{character}, genome is ignored and and empty \code{GRanges} is returned.
#' @param use.cache   (boolean) wether to use or not the cytoband information included in the packge. \code{use.cache=FALSE} will force a download from the UCSC.
#' 
#' @return
#' It returns a \code{\link{GenomicRanges}} object with the cytobands of the specified genome. If no cytobands are available for any reason, an empty \code{GRanges} is returned.
#' 
#' 
#' @note 
#' 
#' This function is memoised (cached) using the \code{\link{memoise}} package. To empty the 
#' cache, use \code{\link{forget}(getCytobands)}
#'  
#' @seealso \code{\link{plotKaryotype}}
#' 
#' @examples
#'  
#' #get the cytobands for hg19 (using the data included in the package)
#' cyto <- getCytobands("hg19")
#' 
#' #do not use the included data and force the download from UCSC
#' cyto <- getCytobands("hg19", use.cache=FALSE)
#' 
#' #get the cytobands for Drosophila Melanogaster
#' cyto <- getCytobands("dm6")
#' 
#' #get the cytobands for Chimpanzee (not included in the package)
#' cyto <- getCytobands("panTro4")
#'  
#' @export getCytobands
#' 


getCytobands <- NULL #Neede so roxygen writes the documentation file


#This is the internal function used. The exported and memoised version is created on-load 
#time as a memoised version of this.

.getCytobands <- function(genome="hg19", use.cache=TRUE) {
  #if it's a custom genome, do not attempt to get the cytobands
  if(!is.character(genome) | length(genome)>1) { 
    return(GRanges())
  }
  
  #If the cytobands are in the included cache, use them
  if(use.cache) {
    if(genome %in% names(data.cache[["cytobands"]])) {
      return(data.cache[["cytobands"]][[genome]])
    }
  }
    
  cytobands <- tryCatch(expr={
    biovizBase::getIdeogram(genome, cytobands=TRUE)
    #Old version. Changed to a dependency on boivizBase as requested by package reviewer.    
    # ucsc.session <- browserSession()
    # genome(ucsc.session) <- genome
    # cytobands <- getTable(rtracklayer::ucscTableQuery(ucsc.session,"cytoBandIdeo"))
    # cytobands$name <- as.character(cytobands$name)
    # cytobands$gieStain <- as.character(cytobands$gieStain)
    # toGRanges(cytobands)    
    },
    error = function(e) {
      message("Error when retrieving from UCSC the Cytobands for ", genome,
              ". Returning no cytobands.", e)
      return(GRanges())
    }
  )
  
  return(cytobands)
  
}

# 
# #Code used to save the predownloaded Cytobands for some common genomes
# cytobands.cache <- list()
# cytobands.cache[["hg19"]] <- getCytobands("hg19")
# cytobands.cache[["hg38"]] <- getCytobands("hg38")
# cytobands.cache[["mm9"]] <- getCytobands("mm9")
# cytobands.cache[["mm10"]] <- getCytobands("mm10")
# cytobands.cache[["rn5"]] <- getCytobands("rn5")
# cytobands.cache[["rn6"]] <- getCytobands("rn6")
# cytobands.cache[["danRer10"]] <- getCytobands("danRer10")
# cytobands.cache[["dm6"]] <- getCytobands("dm6")
# cytobands.cache[["ce6"]] <- GRanges()
# cytobands.cache[["sacCer3"]] <- GRanges()
# 
# genomes.cache <- list()
# genomes.cache[["hg19"]] <- GRangesForUCSCGenome(genome="hg19")
# genomes.cache[["hg38"]] <- GRangesForUCSCGenome(genome="hg38")
# genomes.cache[["mm9"]] <- GRangesForUCSCGenome(genome="mm9")
# genomes.cache[["mm10"]] <- GRangesForUCSCGenome(genome="mm10")
# genomes.cache[["rn5"]] <- GRangesForUCSCGenome(genome="rn5")
# genomes.cache[["rn6"]] <- GRangesForUCSCGenome(genome="rn6")
# genomes.cache[["danRer10"]] <- GRangesForUCSCGenome(genome="danRer10")
# genomes.cache[["dm6"]] <- GRangesForUCSCGenome(genome="dm6")
# genomes.cache[["ce6"]] <- GRangesForUCSCGenome(genome="ce6")
# genomes.cache[["sacCer3"]] <- GRangesForUCSCGenome(genome="sacCer3")
# 
# data.cache <- list(genomes=genomes.cache, cytobands=cytobands.cache)
# devtools::use_data(data.cache, internal = TRUE, overwrite=TRUE)