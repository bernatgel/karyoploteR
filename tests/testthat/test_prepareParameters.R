library(karyoploteR)
context("prepareParameters functions")

test_that("prepareParameters handles missing parameters", {
  
  dd <- data.frame(chr=c("chr1", "chr2", "chr3"),
                   start=c(1, 100, 1000), end=c(1, 200, 10000),
                   y0=c(-1, 0, 1), y1=c(0, 1, 2))
  dd.gr <- toGRanges(dd)  
  
  #create a plot so we get the karyoplot object
  kp <- plotKaryotype()
  
  #Check errors when values are missing. It expects all values to be passed.
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  
  #But most of them can be null or NA
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=NA))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=NA, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=NA, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=NA, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=NA, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=NA, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=NA, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NA, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=NULL))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=NULL, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=0, r1=NULL, ymin=0, ymax=1))
  expect_silent(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=dd$y0, r0=NULL, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=dd$start, y=NULL, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=dd$chr, x=NULL, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=NULL, chr=NULL, x=dd$start, y=dd$y0, r0=0, r1=1, ymin=0, ymax=1))
  
  
  #TODO: Do the same for prepareParameters4
  
})


test_that("prepareParameters has correct argument precedence", {
  
  dd <- data.frame(chr=c("chr1", "chr2", "chr3"),
                   start=c(1, 100, 1000), end=c(1, 200, 10000),
                   y=c(-1, 0, 1), y0=c(-1, 0, 1), y1=c(0, 1, 2))
  dd.gr <- toGRanges(dd)  
  
  #create a plot so we get the karyoplot object
  kp <- plotKaryotype()
  
  #PrepareParameters2
    #chr
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr2", "chr3"))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", "chr1"), r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr1", "chr1"))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", NA), r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr1"))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", NA), r0=0, r1=1, ymin=0, ymax=1, filter.data = FALSE)$chr, c("chr1", "chr1", NA))
    #x
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x, c(1, 150, 5500))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, x=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$x, c(1,2,3))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, x=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$x, c(1,2,NA))
    #y
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y, c(-1, 0, 1))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, y=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$y, c(1,2,3))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, y=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$y, c(1,2,NA))
    
  #PrepareParameters4
    #chr
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr2", "chr3"))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", "chr1"), r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr1", "chr1"))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", NA), r0=0, r1=1, ymin=0, ymax=1)$chr, c("chr1", "chr1"))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr1", NA), r0=0, r1=1, ymin=0, ymax=1, filter.data = FALSE)$chr, c("chr1", "chr1", NA))
    #x0
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x0, c(1, 100, 1000))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x0=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$x0, c(1,2,3))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x0=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$x0, c(1,2,NA))
    #x1
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x1, c(1, 200, 10000))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x1=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$x1, c(1,2,3))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x1=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$x1, c(1,2,NA))
    #y0
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y0, c(-1, 0, 1))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y0=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$y0, c(1,2,3))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y0=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$y0, c(1,2,NA))
    #y1
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y1, c(0,1,2))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y1=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1)$y1, c(1,2,3))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y1=c(1,2,NA), r0=0, r1=1, ymin=0, ymax=1)$y1, c(1,2,NA))
    
})



test_that("prepareParameters returns the correct classes", {
  dd <- data.frame(chr=c("chr1", "chr2", "chr3"),
                   start=c(1, 100, 1000), end=c(1, 200, 10000),
                   y=c(-1, 0, 1), y0=c(-1, 0, 1), y1=c(0, 1, 2))
  dd.gr <- toGRanges(dd)  
  
  #create a plot so we get the karyoplot object
  kp <- plotKaryotype()
  
  #PrepareParameters2
  expect_is(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$chr, "character")
  expect_is(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x, "numeric")
  expect_is(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y, "numeric")
  
  expect_is(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1, filter.data = FALSE)$chr, "character")
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, x=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, y=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  
  #PrepareParameters4
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$chr, "character")
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x0, c("numeric", "integer"))
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y0, c("numeric", "integer"))
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$x1, c("numeric", "integer"))
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y1, c("numeric", "integer"))
  
  expect_is(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c(1,2,3), r0=0, r1=1, ymin=0, ymax=1, filter.data = FALSE)$chr, "character")
  expect_error(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x0=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y0=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, x1=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  expect_error(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, y1=c("1", "1", "1"), r0=0, r1=1, ymin=0, ymax=1))
  
})



test_that("prepareParameters return correct values", {
  
  dd <- data.frame(chr=c("chr1", "chr2", "chr3"),
                   start=c(1, 100, 1000), end=c(1, 200, 10000),
                   y=c(-1, 0, 1), y0=c(-1, 0, 1), y1=c(0, 1, 2))
  dd.gr <- toGRanges(dd)  
  
  #create a plot so we get the karyoplot object
  kp <- plotKaryotype()
  kp2 <- plotKaryotype(chromosomes = "chr1")
  
  #PrepareParameters2
    #chr has already been tested elsewhere and not affected by additional parameters
    #x has already been tested elsewhere and not affected by additional parameters
    #y 
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y, c(-1, 0, 1))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=0.5, ymin=0, ymax=1)$y, c(-0.5, 0, 0.5))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=1, r1=0.5, ymin=0, ymax=1)$y, c(1.5, 1, 0.5))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0.5, r1=0.5, ymin=0, ymax=1)$y, c(0.5, 0.5, 0.5))
    
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=0.5)$y, c(-2, 0, 2))
    expect_equal(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=1, ymax=0.5)$y, c(4, 2, 0))
    expect_true(all(is.infinite(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0.5, ymax=0.5)$y)))

    expect_equal(length(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y), 3)
    expect_equal(length(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2", "chr3", "chr4"), r0=0, r1=1, ymin=0, ymax=1)$y), 4)
    expect_equal(length(prepareParameters2(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2"), r0=0, r1=1, ymin=0, ymax=1)$y), 3)

    expect_equal(length(prepareParameters2(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y), 1)
    expect_equal(length(prepareParameters2(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1, filter.data=FALSE)$y), 3)
    

  #PrepareParameters4
    #chr has already been tested elsewhere and not affected by additional parameters
    #x has already been tested elsewhere and not affected by additional parameters
    #y0
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y0, c(-1, 0, 1))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=0.5, ymin=0, ymax=1)$y0, c(-0.5, 0, 0.5))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=1, r1=0.5, ymin=0, ymax=1)$y0, c(1.5, 1, 0.5))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0.5, r1=0.5, ymin=0, ymax=1)$y0, c(0.5, 0.5, 0.5))
    
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=0.5)$y0, c(-2, 0, 2))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=1, ymax=0.5)$y0, c(4, 2, 0))
    expect_true(all(is.infinite(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0.5, ymax=0.5)$y0)))
    
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y0), 3)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2", "chr3", "chr4"), r0=0, r1=1, ymin=0, ymax=1)$y0), 4)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2"), r0=0, r1=1, ymin=0, ymax=1)$y0), 3)
    
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y0), 1)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1, filter.data=FALSE)$y0), 3)
    
    #y1
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y1, c(0, 1, 2))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=0.5, ymin=0, ymax=1)$y1, c(0, 0.5, 1))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=1, r1=0.5, ymin=0, ymax=1)$y1, c(1, 0.5, 0))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0.5, r1=0.5, ymin=0, ymax=1)$y1, c(0.5, 0.5, 0.5))
    
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=0.5)$y1, c(0, 2, 4))
    expect_equal(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=1, ymax=0.5)$y1, c(2, 0, -2))
    expect_true(all(is.infinite(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0.5, ymax=0.5)$y1)))
    
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y1), 3)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2", "chr3", "chr4"), r0=0, r1=1, ymin=0, ymax=1)$y1), 4)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp, data=dd.gr, chr=c("chr1", "chr2"), r0=0, r1=1, ymin=0, ymax=1)$y1), 3)
    
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1)$y1), 1)
    expect_equal(length(prepareParameters4(function.name = "Test", karyoplot=kp2, data=dd.gr, r0=0, r1=1, ymin=0, ymax=1, filter.data=FALSE)$y1), 3)
    
})    
    
    



