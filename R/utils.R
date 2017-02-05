#internal
#Utility functions used only within the package


#Recycle the arguments as needed.
#Taken from:
# http://stackoverflow.com/questions/9335099/implementation-of-standard-recycling-rules
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

#Colors
#Given a color, returns a lighter one
#TODO: Is there a better way to do that?
lighter <- function(col, amount=150) {
  new.col <- ((grDevices::col2rgb(col))+amount)/255
  new.col[new.col[,1]>1,1] <- 1
  return(grDevices::rgb(t(new.col)))  
}

#Given a color, returns a darker one
darker <- function(col, amount=150) {
  new.col <- ((grDevices::col2rgb(col))-amount)/255
  new.col[new.col[,1]<0, 1] <- 0
  return(grDevices::rgb(t(new.col)))  
}