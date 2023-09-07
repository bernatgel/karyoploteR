#Build Internat Data


#Genomes and Cytobands

# # #Code used to save the predownloaded Cytobands for some common genomes
# genomes <- c("hg18", "hg19", "hg38", "mm9", "mm10", "mm39", "rn5", "rn6", "rn7", "susScr11", 
#              "bosTau9", "bosTau8", "equCab3", "equCab2", "panTro6", "panTro5", "rheMac10",
#              "danRer10", "danRer11", "xenTro10", "dm3", "dm6", 
#              "ce6", "ce10", "ce11", "sacCer2", "sacCer3")
# 
# cytobands.cache <- list()
# genomes.cache <- list()
# 
# for(g in genomes) {
#   message("Downloading data for ", g)
#   cytobands.cache[[g]] <- getCytobands(g, use.cache=FALSE)
#   genomes.cache[[g]] <- GRangesForUCSCGenome(genome=g)
#   Sys.sleep(200)  #To ensure no blocking from UCSC due to too many requests
# }
# 
# data.cache <- list(genomes=genomes.cache, cytobands=cytobands.cache)
# 
# library(devtools)
# use_data(data.cache, internal = TRUE, overwrite=TRUE)
#  
# load("R/sysdata.rda")
# data.cache

#OLD CODE
# cytobands.cache <- list()
# cytobands.cache[["hg19"]] <- getCytobands("hg19", use.cache=FALSE)
# cytobands.cache[["hg38"]] <- getCytobands("hg38", use.cache=FALSE)
# cytobands.cache[["mm9"]] <- getCytobands("mm9", use.cache=FALSE)
# cytobands.cache[["mm10"]] <- getCytobands("mm10", use.cache=FALSE)
# cytobands.cache[["rn5"]] <- getCytobands("rn5", use.cache=FALSE)
# cytobands.cache[["rn6"]] <- getCytobands("rn6", use.cache=FALSE)
# cytobands.cache[["danRer10"]] <- getCytobands("danRer10", use.cache=FALSE)
# cytobands.cache[["dm6"]] <- getCytobands("dm6", use.cache=FALSE)
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

#Character Polygons for kpChars and kpPlotSequence
#Download https://github.com/heike/gglogo/raw/master/data/alphabet.rda
#load("alphabet.rda")
#char.polygons <- split(alphabet, alphabet$group)

# library(devtools)
# use_data(data.cache, char.polygons, internal = TRUE, overwrite=TRUE)
# 
# load("R/sysdata.rda")
# data.cache
