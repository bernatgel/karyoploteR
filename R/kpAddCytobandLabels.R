#' kpAddCytobandLabels
#' 
#' @description 
#' 
#' Plots the base numbers along the chromosome ideograms
#' 
#' @details 
#'  
#' This function can be used to add labels idenfifying the cytobands. It gets the
#' labels from the cytobands information stored in the karyoplot object and it will
#' only plot the labels that fit inside the available space. This means than in some 
#' cases (such as when plotting a complete genome with default parameters) it is 
#' possible that no labels at all are added.
#' 
#' @usage kpAddCytobandLabels(karyoplot, cex=0.5, force.all=FALSE, ...)
#' 
#' @param karyoplot  (karyoplot object) A valid karyoplot object created by a call to \code{\link{plotKaryotype}}
#' @param cex  (numeric) The cex parameter for the cytoband labels
#' @param force.all  (boolean) If true, all cytoband labels are plotted, even if they do not fit into the cytobands (Defaults to FALSE)
#' @param ...  Any other parameter to be passed to internal function calls. Specially useful for graphic parameters.
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{plotKaryotype}}
#' 
#' @examples
#'   
#' kp <- plotKaryotype()
#' kpAddBaseNumbers(kp)
#' kpAddCytobandLabels(kp)
#' 
#' kp <- plotKaryotype(chromosomes="chr17")
#' kpAddBaseNumbers(kp, tick.dist=10000000, minor.tick.dist=1000000)
#' kpAddCytobandLabels(kp)
#' 
#'  
#' @export kpAddCytobandLabels
#' 

kpAddCytobandLabels <- function(karyoplot, cex=0.5, force.all=FALSE,  ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
      
  if(!is.null(karyoplot$cytobands)) {
    if(length(karyoplot$cytobands)>0) { #If there are cytobands to plot, plot them
  
      labels <- karyoplot$cytobands$name
      bands.chr <- as.character(seqnames(karyoplot$cytobands))
      
      ylabel <- karyoplot$ideogram.mid(bands.chr)
      
      bandxleft <- ccf(x=start(karyoplot$cytobands), chr=bands.chr)$x
      bandxright <- ccf(x=end(karyoplot$cytobands), chr=bands.chr)$x
      
      bandmids <- (bandxleft +(bandxright-bandxleft)/2)
      
      label.width <- graphics::strwidth(labels, cex=cex)
      
      if(force.all == FALSE) {
        do.fit <- label.width < (bandxright - bandxleft)
      } else {
        do.fit <- rep(TRUE, length(labels))
      }    
      
      if(any(do.fit)) {
        graphics::text(x=bandmids[do.fit], y=ylabel[do.fit], labels=labels[do.fit], cex=cex, ...)
      }
      
    }
  }
  #If no cytobands are available, do nothing
  
  invisible(karyoplot)
}
