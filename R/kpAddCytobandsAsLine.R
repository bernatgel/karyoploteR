#' kpAddCytobandsAsLine
#' 
#' @description 
#' 
#' Plots the chromosome cytobands in a karyoplot as a line
#' 
#' @details 
#' 
#' Plots the cytobands representing the chromosome structure in a karyoplot. It extracts the 
#' cytobands from the \code{karyoplot} object it recieves as a parameter. It is possible to 
#' specify the colors used to plot the cytobands. In contrast to \code{\link{kpAddCytobands}}
#' it represents the chromosomes as a thin line
#' 
#' @note In general, this function is automatically called by plotKaryotype
#' and the user never nees to call it. 
#' 
#' @usage kpAddCytobandsAsLine(karyoplot, color.table=NULL, color.schema='only.centromeres', lwd=3, lend=1, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param color.table  (named character vector) a table specifying the colors to plot the cytobands. If NULL, it gets the colors calling \code{getCytobandColors}. (defaults to NULL)
#' @param color.schema  (character: 'only.centromeres', 'circos', 'biovizbase') The name of the color schema to use. It is directly passed along to \code{\link{getCytobandColors}}. \code{color.table} takes precendence over \code{color.schema}. (defaults to 'only.centromeres')
#' @param lwd (integer) The width of the line used to represent the ideogram (defaults to 3)
#' @param lend (0, 1 or 2) The type of line end. (defaults to 1, "butt")
#' @param ...  any additional parameter to be passed to the functions called from kpAddCytobands.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{getCytobandColors}}, \code{\link{kpAddBaseNumbers}}, \code{\link{kpAddCytobandLabels}}
#' 
#' @examples
#'
#' 
#' kp <- plotKaryotype(ideogram.plotter = NULL)
#' kpAddCytobandsAsLine(kp)
#'  
#' @export kpAddCytobandsAsLine
#' 



kpAddCytobandsAsLine <- function(karyoplot, color.table=NULL, color.schema='only.centromeres', lwd=3, lend=1, ...) {
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  color.table <- getCytobandColors(color.table, color.schema)
 
  
  #If there are no cytobands, create a single "fake" cytoband to represent the whole chromosome
  if(!is.null(karyoplot$cytobands) && length(karyoplot$cytobands)>0) { 
    cyto <- karyoplot$cytobands  
  } else {
    cyto <- karyoplot$genome
    mcols(cyto) <- data.frame(name=seqnames(cyto), gieStain="gpos50", stringsAsFactors=FALSE)  
  }
  
  #filter out the cytobands out of our chromosomes
  cyto <- filterChromosomes(cyto, keep.chr = karyoplot$chromosomes)
  
  #And plot them
  ybottom <- mids(as.character(seqnames(cyto)))
  ytop <- mids(as.character(seqnames(cyto)))
    
  xleft <- ccf(x=start(cyto), chr=as.character(seqnames(cyto)))$x
  xright <- ccf(x=end(cyto), chr=as.character(seqnames(cyto)))$x
    
  col <- color.table[as.character(cyto$gieStain)]
    
  graphics::segments(x0 = xleft, x1=xright, y0=ybottom, y1=ytop, col=col, lwd=lwd, lend=lend)      

  
  invisible(karyoplot)
}

