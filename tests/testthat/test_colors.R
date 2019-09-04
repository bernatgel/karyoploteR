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
  

test_that("darker, lighter and transparent work correctly", {
  
  expect_equal(darker(c("red", "#00FF00", NA, "black")), c("#690000", "#006900", NA, "#000000"))
  expect_equal(darker(c("red", "#00FF00", NA, "black"), amount=20), c("#EB0000", "#00EB00", NA, "#000000"))
  expect_error(darker(c("red", "#00FF00", NA, "black"), amount=-1), "amount must be a number between 0 and 255")
  expect_error(darker(c("red", "#00FF00", NA, "black"), amount=FALSE), "amount must be a single number")
  expect_error(darker(c("red", "#00FF00", NA, "black"), amount="GREEN"), "amount must be a single number")
  
  expect_equal(lighter(c("red", "#00FF00", NA, "white")), c("#FF9696", "#96FF96", NA, "#FFFFFF"))
  expect_equal(lighter(c("red", "#00FF00", NA, "white"), amount=20), c("#FF1414", "#14FF14", NA, "#FFFFFF"))
  expect_error(lighter(c("red", "#00FF00", NA, "black"), amount=-1), "amount must be a number between 0 and 255")
  expect_error(lighter(c("red", "#00FF00", NA, "black"), amount=FALSE), "amount must be a single number")
  expect_error(lighter(c("red", "#00FF00", NA, "black"), amount="GREEN"), "amount must be a single number")
  
  expect_equal(transparent(c("red", "#00FF00", NA, "black")), c("#FF000080", "#00FF0080", NA, "#00000080"))
  expect_equal(transparent(c("red", "#00FF00", NA, "white"), amount=0.2), c("#FF0000CC", "#00FF00CC", NA, "#FFFFFFCC"))
  expect_error(transparent(c("red", "#00FF00", NA, "black"), amount=-1), "amount must be a number between 0 and 1")
  expect_error(transparent(c("red", "#00FF00", NA, "black"), amount=FALSE), "amount must be a single number")
  expect_error(transparent(c("red", "#00FF00", NA, "black"), amount="GREEN"), "amount must be a single number")
   
})  
    



