#' plotCytobands
#' 
#' @description 
#' 
#' Plots the chromosome cytobands in a karyoplot
#' 
#' @details 
#' 
#' Plots the cytobands representing the chromosome structure in a karyoplot. It extracts the 
#' cytobands from the \code{karyoplot} object it recieves as a parameter. It is possible to 
#' specify the colors used to plot the cytobands. 
#' 
#' @note In general, this function is automatically called by plotKaryotype
#' and the user never nees to call it. 
#' 
## @usage plotCytobands(karyoplot, color.table=NULL, add.cytobands.names=FALSE, add.base.numbers=FALSE, ...)
#' @usage plotCytobands(karyoplot, color.table=NULL, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param color.table  (named character vector) a table specifying the colors to plot the cytobands. If NULL, it gets the colors calling \code{getCytobandColors}. (defaults to NULL)
## @param add.cytobands.names  (boolean) whether to add or not the cytoband names to the plot. (defaults to FALSE)
## @param add.base.numbers  (boolean) whether to add the base numbers to the plot. (defaults to FALSE)
#' @param ...  any additional parameter to be passed to the functions called from plotCytobands.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{getCytobandColors}}, \code{\link{kpAddBAseNumbers}}, \code{\link{kpAddCytobandNames}}
#' 
#' @examples
#'
#' 
#' kp <- plotKaryotype(ideogram.plotter = NULL)
#' plotCytobands(kp)
#'  
#' @export plotCytobands
#' 



plotCytobands <- function(karyoplot, color.table=NULL, add.cytobands.names=FALSE,
                          add.base.numbers=FALSE, ...) {
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
    
  color.table <- getCytobandColors(color.table)
  
  #If there are no cytobands, create a single "fake" cytoband to represent the whole chromosome
  if(!is.null(karyoplot$cytobands) && length(karyoplot$cytobands)>0) { 
    cyto <- karyoplot$cytobands  
  } else {
    cyto <- karyoplot$genome
    mcols(cyto) <- data.frame(name=seqnames(cyto), gieStain="gpos50", stringsAsFactors=FALSE)  
  }
  
  #And plot them  
  
  ybottom <- mids(as.character(seqnames(cyto))) - pp$ideogramheight/2
  ytop <- mids(as.character(seqnames(cyto))) + pp$ideogramheight/2
    
  xleft <- ccf(x=start(cyto))$x
  xright <- ccf(x=end(cyto))$x
    
  #col <- do.call(c, color.table[cyto$gieStain])
  col <- color.table[as.character(cyto$gieStain)]
    
  graphics::rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col)      
  #   
  # if(add.cytobands.names) {
  #   plotCytobandsLabels(karyoplot=karyoplot, ...)
  # }
  #   
  # if(add.base.numbers) {
  #   kpAddBaseNumbers(karyoplot, ...)
  # }
  
  invisible(karyoplot)
}

