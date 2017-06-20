library(karyoploteR)
context("Cytobands and Cytoband Labels")

#Test getCytobands
test_that("getCytobands returns the expected cytobands", {
  cytobands.hg19 <- getCytobands("hg19")
  expect_is(cytobands.hg19, "GRanges")
  expect_equal(names(mcols(cytobands.hg19)), c("name", "gieStain"))
  
  cytobands.hg19 <- filterChromosomes(cytobands.hg19, chr.type = "canonical", organism = "hg")
  expect_equal(length(cytobands.hg19), 862)
  
  cytobands.hg19.chr17 <- filterChromosomes(cytobands.hg19, keep.chr = "chr17")
  expect_equal(cytobands.hg19.chr17$name, 
               c('p13.3','p13.2','p13.1','p12','p11.2','p11.1','q11.1','q11.2','q12','q21.1','q21.2','q21.31','q21.32','q21.33','q22','q23.1','q23.2','q23.3','q24.1','q24.2','q24.3','q25.1','q25.2','q25.3')
  )

  #Test with other valid genomes
  expect_equal(length(getCytobands("mm10")), 448)
  expect_equal(length(getCytobands("rn6")), 953)
  expect_equal(length(getCytobands("dm6")), 6917)
  
  
  #Check null values  
  expect_equal(length(getCytobands(NULL)), 0)
  expect_silent(getCytobands(NULL))
  expect_equal(length(getCytobands(genome=NA)), 0)
  expect_silent(getCytobands(genome=NA))
  expect_equal(length(getCytobands(genome=NULL)), 0)
  expect_silent(getCytobands(genome=NULL))
  
  #Check with invalid genomes
  #FIXIT
  #Should not fail with error, but message and return empty object:  expect_message(getCytobands(genome="INVALID_GENOME"))
  
  
})

#kpAddCytobands
test_that("kpAddCytobands does not error. Not checking correct plotting.", {
  #Note: We are not testing that the generated image is correct!
  
  #all standard
  kp <- plotKaryotype()
  expect_silent(kpAddCytobandLabels(kp))
  expect_silent(kpAddCytobandLabels(kp, force.all = TRUE, cex = 1, srt=90))
    
  #with only some chromosomes
  kp <- plotKaryotype(chromosomes = "chr2")
  #FIXIT: expect_silent(kpAddCytobandLabels(kp))
  expect_silent(kpAddCytobandLabels(kp, force.all = TRUE, cex = 1, srt=90))
  
  
  #with a custom genome
  custom.genome <- toGRanges(data.frame(chr="A", start=0, end=1000))
  kp <- plotKaryotype(genome=custom.genome)
  expect_silent(kpAddCytobandLabels(kp))
  
  
})


