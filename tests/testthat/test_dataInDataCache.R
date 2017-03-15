library(karyoploteR)
context("Content of the data cache")

#Test order of genomes
test_that("the genomes are sorted", {
  dc <- karyoploteR:::data.cache
  
  #hg19, canonical chromosomes
  expect_equal(paste0(seqlevels(dc[["genomes"]][["hg19"]])[1:24], collapse = ","), paste0(paste0("chr", c(1:22, "X", "Y")), collapse = ","))
  #hg38, canonical chromosomes
  expect_equal(paste0(seqlevels(dc[["genomes"]][["hg19"]])[1:24], collapse = ","), paste0(paste0("chr", c(1:22, "X", "Y")), collapse = ","))  

  #TODO: Check the other genomes
  
})
