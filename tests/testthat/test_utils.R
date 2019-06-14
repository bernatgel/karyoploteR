library(karyoploteR)
context("utils.R functions")


test_that("preprocess_r0_r1 works as expected", {
  
  #Create the needed kp
  kp <- plotKaryotype()
  
  #Correct specifications
  expect_equal(preprocess_r0_r1(kp, r0=0.2, r1=0.8, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=list(r0=0.2, r1=0.8), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=c(r0=0.2, r1=0.8), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=c(r1=0.8, r0=0.2), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=c(Other.name=2, r1=0.8, r0=0.2), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=c(0.2, 0.8), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=c(0.2, 0.8, r4=4), r1=NULL, data.panel = 1), list(r0=0.2, r1=0.8))
  
  #Correct including NULL
  expect_equal(preprocess_r0_r1(kp, r0=NULL, r1=0.8, data.panel = 1), list(r0=0, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=NULL, r1=NULL, data.panel = 1), list(r0=0, r1=1))
  expect_equal(preprocess_r0_r1(kp, r0=0.2, r1=NULL, data.panel = 1), list(r0=0.2, r1=1))
  expect_equal(preprocess_r0_r1(kp, r0=list(r0=NULL, r1=0.8), r1=NULL, data.panel = 1), list(r0=0, r1=0.8))
  expect_equal(preprocess_r0_r1(kp, r0=list(r0=NULL, r1=NULL), r1=NULL, data.panel = 1), list(r0=0, r1=1))
  expect_equal(preprocess_r0_r1(kp, r0=list(r0=0.2, r1=NULL), r1=NULL, data.panel = 1), list(r0=0.2, r1=1))
  
  
  #Errors
  #TODO: some of these will error on condition (for is.na) longer than 1 in if. Could return better error messages
  expect_error(preprocess_r0_r1(kp, r0=NA, r1=0.8, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=NA, r1=NA, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=0.2, r1=NA, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=list(r0=NA, r1=0.8), r1=NA, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=list(r0=NA, r1=NA), r1=NA, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=list(r0=0.2, r1=NA), r1=NA, data.panel = 1))
  
  expect_error(preprocess_r0_r1(kp, r0=list(r0=0.2, r1=0.8), r1=1, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=c(r0=0.2, r1=0.8), r1=1, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=c(r0="a", r1=0.8), r1=1, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=list(r0="a", r1=0.8), r1=1, data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=c(r0=0.2, r1=0.8), r1="a", data.panel = 1))
  expect_error(preprocess_r0_r1(kp, r0=list(r0=0.2, r1=0.8), r1="a", data.panel = 1))
  
    
})    
    
    



