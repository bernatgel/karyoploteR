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




#' filterParams
#' 
#' @description 
#' Given a list, select just only the valid.elements from each member. Also 
#' works with vectors instead of lists
#' 
#' @details 
#' This function is used in filtering the graphical parameters when plotting
#' only a part of the genome. For each element of the list, if it has the 
#' exact specified length, filters it using the 'valid.elements' parameter.
#' 
#' @usage filterParams(p, valid.elements, orig.length)
#' 
#' @param p a list or a single vector
#' @param valid.elements a boolean vector with the elements to keep
#' @param orig.length the length of the elements on which to apply the filtering
#' 
#' @return
#' p with some members filtered
#' 
#' 
#' @examples
#'  
#' a <- 1:10
#' b <- 3:5
#' c <- 2
#' 
#' filterParams(list(a,b,c), c(rep(TRUE,5), rep(FALSE,5)), 10)
#' filterParams(a, c(rep(TRUE,5), rep(FALSE,5)), 10)
#'  
#' @export filterParams
#' 
filterParams <- function(p, valid.elements, orig.length) {
  if(methods::is(p, "list")) { #If p is a list, filter each element independently
    for(i in seq_len(length(p))) {
      if(length(p[[i]])==orig.length) {
        p[[i]] <- p[[i]][valid.elements]
      }
    }
  } else { #else, filter p as a single element
    if(length(p)==orig.length) {
      p <- p[valid.elements]
    }
  }
  return(p)
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
#' @param col (color) The original color
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
#' @param col (color) The original color
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



#Autotrack
#TODO: document and export? maybe export only if used out of prepareParameters2 
#      and 4

processAutotrack <- function(r0, r1, autotrack) {
  if(!all(unlist(lapply(autotrack, methods::is, "numeric")))) stop("'autotrack' must be a list of numerics and it is a ", unlist(lapply(autotrack, class)))
  if(length(autotrack)<2) stop("'autotrack' must be a list of numerics of length 2 or 3 and it is of length ", length(autotrack))    
  
  at.current.min <- min(autotrack[[1]])
  at.current.max <- max(autotrack[[1]])
  at.total <- autotrack[[2]]
  at.margin <- ifelse(length(autotrack)==2, 0.05, autotrack[[3]])
  
  tr.height <- (r1-r0)/at.total
  r0 <- r0+(at.current.min-1)*tr.height
  r1 <- r0+(at.current.max-at.current.min+1)*tr.height-tr.height*at.margin
  message("r0=", r0, "   r1=", r1)
  return(c(r0=r0, r1=r1))
}