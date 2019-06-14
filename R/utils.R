#internal
#Utility functions used only within the package


#Recycle the arguments as needed.
#Taken from:
# http://stackoverflow.com/questions/9335099/implementation-of-standard-recycling-rules
recycle <- function(...){
  dotList <- list(...)
  max.length <- max(vapply(dotList, length, FUN.VALUE=0))
  lapply(dotList, rep, length=max.length)
}


#Only recycles the first argument and returns it
recycle.first <- function(...){
  dotList <- list(...)
  max.length <- max(vapply(dotList, length, FUN.VALUE = 0))
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


############  Clipping   ###############
#' processClipping
#' 
#' @description 
#' Sets image clipping if needed
#' 
#' @details 
#' Small utility function to help manage clipping. If the current plot is
#' a zoomed plot and clipping is TRUE, activate the clip to the current 
#' data.panel. This will hide any plotting ocurring out of the data.panel
#' region.
#' 
#' @note Users wont usually use this function. It is used by the plotting functions 
#' to set the clipping if needed
#' 
#' @usage processClipping(karyoplot, clipping, data.panel) 
#' 
#' @param karyoplot (KaryoPlot) A KaryoPlot object representing the current plot
#' @param clipping (logical) Wheter clipping should be activated or not
#' @param data.panel (data panel identifier) The name of the data panel on which the plot should be allowed. Anything plotted outside it will be hidden (if clipping==TRUE and the plot is a zoom plot)
#' 
#' @return
#' Returns the original karyoplot object, unchanged.
#' 
#' @examples
#' 
#' kp <- plotKaryotype()
#' processClipping(kp, TRUE, 1)
#'
#'  
#' @export processClipping
#'
processClipping <- function(karyoplot, clipping, data.panel) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")

  if(karyoplot$zoom==TRUE) {
    if(clipping==TRUE) {
      dpbb <- getDataPanelBoundingBox(karyoplot, data.panel)
      graphics::clip(x1 = dpbb$x0, x2 = dpbb$x1, y1 = dpbb$y0, y2=dpbb$y1)
    }
  }
  invisible(karyoplot)
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





#Internal function. Not exported. Used to validate and preprocess r0 and r1
preprocess_r0_r1 <- function(karyoplot, r0, r1, data.panel) {
  
  if(!is.null(r0) && is.null(r1)) { #Maybe r0 contains the r0 and r1 information
    #It might be a list
    if(is.list(r0) && 
       utils::hasName(r0, "r0") && utils::hasName(r0, "r1") &&
       (is.null(r0$r0) || is.numeric(r0$r0)) && (is.null(r0$r1) || is.numeric(r0$r1))) {
      r1 <- r0$r1
      r0 <- r0$r0
    } else {
      #It might be a two element array
      if(is.numeric(r0) && length(r0)>=2) {
        if(all(c("r0", "r1") %in% names(r0))) {
          r1 <- setNames(r0["r1"], NULL)
          r0 <- setNames(r0["r0"], NULL)
        } else {
          r1 <- setNames(r0[2], NULL)
          r0 <- setNames(r0[1], NULL)
        }
      }
    }
  } 
  if(is.null(r0)) {
    r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  }
  if(is.null(r1)) {
    r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  }
  
  
  #Finally, check that r0 and r1 are valid numbers
  if(!(is.numeric(r0) && length(r0)==1)) stop("Invalid r0 specification. Check karyoploteR's documentation for more information.")
  if(!(is.numeric(r1) && length(r1)==1)) stop("Invalid r1 specification. Check karyoploteR's documentation for more information.")
  
  return(list(r0=r0, r1=r1))  
}
