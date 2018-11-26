library(karyoploteR)
context("Main plotKaryotype function")

#Test plotKaryotype
test_that("the function fails as expected with invalid parameters", {
  
  valid.genome <- toGRanges(data.frame(chr=c("1"), start=c(1), end=1000))
  invalid.genome <- toGRanges(data.frame(chr=c("1"), start=c(1,1), end=1000))
  
  expect_silent(plotKaryotype())
  expect_silent(plotKaryotype(genome=valid.genome))
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
