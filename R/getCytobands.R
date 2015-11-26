



getCytobands <- memoise(function(genome="hg19", force.download=FALSE) {
  if(!is.character(genome) | length(genome)>1) { #if it's a custom genome, do not attempt to get the cytobands
    return(GRanges())
  }
  
  if(!force.download) { #If the cytobands are in the included cache, use them
    if(genome %in% names(data.cache[["cytobands"]])) {
      return(data.cache[["cytobands"]][[genome]])
    }
  }
    
  cytobands <- tryCatch(expr={
    ucsc.session <- browserSession()
    genome(ucsc.session) <- genome
    cytobands <- getTable(ucscTableQuery(ucsc.session,"cytoBandIdeo"))
    cytobands$name <- as.character(cytobands$name)
    cytobands$gieStain <- as.character(cytobands$gieStain)
    toGRanges(cytobands)    
  },
    error = function(e) {
      message(paste0("Error when retrieving from UCSC the Cytobands for ", genome,". Returning no cytobands.", e))
      return(GRanges())
    })
  
  return(cytobands)
  
})

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