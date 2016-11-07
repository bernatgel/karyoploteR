#internal
#Utility functions used only within the package


#Recycle the arguments as needed.
#Taken from: http://stackoverflow.com/questions/9335099/implementation-of-standard-recycling-rules
recycle <- function(...){
  dotList <- list(...)
  max.length <- max(sapply(dotList, length))
  lapply(dotList, rep, length=max.length)
}


#Only recycles the first argument and returns it
recycle.first <- function(...){
  dotList <- list(...)
  max.length <- max(sapply(dotList, length))
  return(rep_len(dotList[[1]], length.out=max.length))
}
