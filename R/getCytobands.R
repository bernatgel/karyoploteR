



getCytobands <- memoise(function(genome="hg19") {
  if(!is.character(genome) | length(genome)>1) { #if it's a custom genome, do not attempt to get the cytobands
    return(GRanges())
  }
  
  ucsc.session <- browserSession()
  genome(ucsc.session) <- genome
  cytobands <- getTable(ucscTableQuery(ucsc.session,"cytoBandIdeo"))
    
  cytobands$name <- as.character(cytobands$name)
  cytobands$gieStain <- as.character(cytobands$gieStain)
  
  
  return(toGRanges(cytobands))
  
})

