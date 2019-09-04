library(karyoploteR)
context("colors.R functions")


test_that("preprocessColors works correctly", {
  
  expect_equal(preprocessColors(col="red", border="green"), list(col="red", border="green"))
  expect_equal(preprocessColors(col="red", border=NA), list(col="red", border=NA))
  expect_equal(preprocessColors(col=NA, border="green"), list(col=NA, border="green"))
  expect_equal(preprocessColors(col=NA, border=NA), list(col=NA, border=NA))
  expect_equal(preprocessColors(col=c(NA, "red"), border=c("green", NA)), list(col=c(NA, "red"), border=c("green", NA)))
  expect_equal(preprocessColors(col=c(NA, "red"), border=NULL), list(col=c(NA, "red"), border=c(NA, "#9B0000")))
  expect_equal(preprocessColors(col=NULL, border=c("green", NA)), list(col=c("#64FF64", NA), border=c("green", NA)))
  expect_equal(preprocessColors(col=NULL, border=NULL), list(col="gray70", border="#4F4F4F"))
  
  
})    
    
    



