



getCytobands <- memoise(function(genome="hg19") {
  if(!is.character(genome) | length(genome)>1) {
    return(GRanges())
  }
  
  
  ucsc.session <- browserSession()
  genome(ucsc.session) <- genome
  cytobands <- getTable(ucscTableQuery(session,"cytoBandIdeo"))
    
  return(toGRanges(cytobands))
  
})

forget(f=getCytobands)