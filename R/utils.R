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

############  Colors  ###############
#' lighter
#' 
#' @description 
#' Given a color, return a lighter one
#' 
#' @details 
#' Very simple utility function to create lighter colors. Given a color, it
#' transforms it to rgb space, adds a set amount to all chanels and transforms
#' it back to a color.
#' 
#' @usage lighter(col, amount=150)
#' 
#' @param color (color) The original color
#' @param amount (integer, [0-255]) The fixed amount to add to each RGB channel (Defaults to 150).
#' 
#' @return
#' A lighter color
#' 
#' @seealso \code{\link{darker}}
#' 
#' @examples
#'  
#' lighter("red")
#' lighter("#333333")
#'  
#' @export lighter
#' 

lighter <- function(col, amount=150) {
  new.col <- ((grDevices::col2rgb(col))+amount)/255
  new.col[new.col[,1]>1,1] <- 1
  return(grDevices::rgb(t(new.col)))  
}

#' darker
#' 
#' @description 
#' Given a color, return a darker one
#' 
#' @details 
#' Very simple utility function to create darker colors. Given a color, it
#' transforms it to rgb space, adds a set amount to all chanels and transforms
#' it back to a color.
#' 
#' @usage darker(col, amount=150)
#' 
#' @param color (color) The original color
#' @param amount (integer, [0-255]) The fixed amount to add to each RGB channel (Defaults to 150).
#' 
#' @return
#' A darker color
#' 
#' @seealso \code{\link{lighter}}
#' 
#' @examples
#'  
#' darker("red")
#' darker("#333333")
#'  
#' @export darker
#'

#Given a color, returns a darker one
darker <- function(col, amount=150) {
  new.col <- ((grDevices::col2rgb(col))-amount)/255
  new.col[new.col[,1]<0, 1] <- 0
  return(grDevices::rgb(t(new.col)))  
}