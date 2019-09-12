library(karyoploteR)
context("Main plotKaryotype function")

#Test plotKaryotype
test_that("plotKaryotype works with different genome definitions", {
  
  valid.genome <- toGRanges(data.frame(chr=c("1"), start=c(1), end=1000))
  invalid.genome <- toGRanges(data.frame(chr=c("1"), start=c(1,1), end=1000)) #repeted chr names
  
  require(BSgenome.Hsapiens.UCSC.hg19)
  
  expect_silent(plotKaryotype())
  expect_silent(plotKaryotype(genome="hg19"))
  expect_silent(plotKaryotype(genome="mm10"))
  expect_silent(plotKaryotype(genome=BSgenome.Hsapiens.UCSC.hg19))
  expect_silent(plotKaryotype(genome=seqinfo(BSgenome.Hsapiens.UCSC.hg19)))
  expect_silent(plotKaryotype(genome=valid.genome))
  expect_silent(plotKaryotype(genome=valid.genome, cex=2))
  expect_silent(plotKaryotype(genome=valid.genome, zoom=toGRanges(data.frame("1", 5, 20))))
  expect_silent(plotKaryotype(genome=valid.genome, ideogram.plotter = NULL))
  expect_silent(plotKaryotype(genome=valid.genome, labels.plotter = NULL))
  expect_error(plotKaryotype(genome=invalid.genome))
  
})

#kpAddCytobands
test_that("the returned karyoplot object is complete and correct", {
  #TODO!
  
})


#kpAddBaseNumbers
test_that("kpAddBaseNumbers works", {
  #github issue 26
  kp <- plotKaryotype(genome="hg19",plot.type=2)
  kpAddBaseNumbers(kp)
  #TODO: Extend testing
  
})


#kpAddCytobandLabels
test_that("kpAddCytobandLabels works", {
  #identified when fixing github issue 26
  kp <- plotKaryotype(genome="hg19", plot.type=2)
  kpAddCytobandLabels(kp)
  
  #TODO: Extend testing
  
})
