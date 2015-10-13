



getCytobands <- memoise(function(genome="hg19") {
  if(!is.character(genome) | length(genome)>1) {
    return(GRanges())
  }
  
  
  ucsc.session <- browserSession()
  genome(ucsc.session) <- genome
  cytobands <- getTable(ucscTableQuery(session,"cytoBandIdeo"))
    
  cytobands$name <- as.character(cytobands$name)
  cytobands$gieStain <- as.character(cytobands$gieStain)
  
  
  return(toGRanges(cytobands))
  
})

forget(f=getCytobands)