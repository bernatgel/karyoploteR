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






############  Autotrack  ###############

#' autotrack
#' 
#' @description 
#' Computes r0 and r1 given track definition
#' 
#' @details 
#' Small utility function to help compute r0 and r1 given the total number of tracks 
#' and the track(s) the current plot will occupy. It also takes into account a margin
#' between tracks and original r0 and r1, so we can say something like, "Out of 5 
#' tracks between 0 and 0.5, this plot will be at track 2", and it will return
#' r0=0.1 and r1=0.2
#' 
#' 
#' @usage autotrack(current.track, total.tracks, margin=0.05, r0=0, r1=1)
#' 
#' @param current.track (numeric) The track or tracks the current plot will occupy, starting from 1. If more than one value is provided, the plot will expand from min(current.track) to max(current.track).
#' @param total.tracks (numeric) The total number of tracks
#' @param margin (numeric) The margin is specified as the part of a track, by default 0.05, 5 percent of the track height. 
#' @param r0 (numeric) the original r0
#' @param r1 (numeric) the original r1
#' 
#' @return
#' A list of two numerics: r0 and r1
#' 
#' @examples
#'
#' #first track out of 4  
#' autotrack(1, 4)
#'
#' #the same, but without margin
#' autotrack(1, 4, 0)
#' 
#' #first and second tracks out of 4
#' autotrack(c(1,2), 4)
#' 
#' #The first track out of 4, fitting the four track between 0 and 0.5
#' autotrack(1, 4, r0=0, r1=0.5)
#'  
#' @export autotrack
#'

autotrack <- function(current.track, total.tracks, margin=0.05, r0=0, r1=1) {
  if(!methods::is(current.track, "numeric")) stop("current.track must be numeric")
  if(!methods::is(total.tracks, "numeric")) stop("total.tracks must be numeric")
  if(length(total.tracks)>1) {
    warning("total.tracks has more than one value. Using only the first one.")
    total.tracks <- total.tracks[1]  
  }
  if(!methods::is(margin, "numeric")) stop("margin must be numeric")
  if(length(margin)>1) {
    warning("margin has more than one value. Using only the first one.")
    margin <- margin[1]  
  }
  if(!methods::is(r0, "numeric")) stop("r0 must be numeric")
  if(!methods::is(r1, "numeric")) stop("r1 must be numeric")
  
  
  at.current.min <- min(current.track)
  at.current.max <- max(current.track)
  
  tr.height <- (r1-r0)/total.tracks
  r0 <- r0+(at.current.min-1)*tr.height
  r1 <- r0+(at.current.max-at.current.min+1)*tr.height-tr.height*margin
  
  return(list(r0=r0, r1=r1))
}


