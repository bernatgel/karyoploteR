
#Note: The function getCytobands is a memoised function. It's better to define it on package
#load so it does not get linked to the memoise version available at installation time.
#See comment by @mtmorgan at https://github.com/Bioconductor/Contributions/issues/199

getCytobands <- NULL

.onLoad <- function(libname, pkgname) {
  getCytobands <<- memoise::memoise(.getCytobands)
}