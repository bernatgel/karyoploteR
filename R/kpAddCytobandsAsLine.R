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
#' and the user never needs to call it. 
#' 
#' @usage kpAddCytobandsAsLine(karyoplot, color.table=NULL, color.schema='only.centromeres', lwd=3, lend=1, clipping=TRUE, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param color.table  (named character vector) a table specifying the colors to plot the cytobands. If NULL, it gets the colors calling \code{getCytobandColors}. (defaults to NULL)
#' @param color.schema  (character: 'only.centromeres', 'circos', 'biovizbase') The name of the color schema to use. It is directly passed along to \code{\link{getCytobandColors}}. \code{color.table} takes precendence over \code{color.schema}. (defaults to 'only.centromeres')
#' @param lwd (integer) The width of the line used to represent the ideogram (defaults to 3)
#' @param lend (0, 1 or 2) The type of line end. (defaults to 1, "butt")
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, cytoband representation will be not drawn out of the drawing are (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the cytobands representation may overflow into the margins of the plot. (defaults to TRUE)
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

#TODO: Since most of the code (and growing) is shared between this function and 
#kpAddCytobands, merge them in an internal function and keep the two separate 
#functions in the public API.

#TODO: Add a cytobands parameter to allow explicitly passing the cytobands 
#in the function call instead of relying in they being available in the 
#KaryoPlot object

kpAddCytobandsAsLine <- function(karyoplot, color.table=NULL, color.schema='only.centromeres', lwd=3, lend=1, clipping=TRUE, ...) {
  
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  #If there are no cytobands, create a single "fake" cytoband to represent the whole chromosome
  if(!is.null(karyoplot$cytobands) && length(karyoplot$cytobands)>0) { 
    cyto <- karyoplot$cytobands  
  } else {
    cyto <- karyoplot$genome
    mcols(cyto) <- data.frame(name=seqnames(cyto), gieStain="gpos50", stringsAsFactors=FALSE)  
  }
  
  if(!methods::is(cyto, "GRanges")) stop("'cytobands' must be a GRanges object")

  
  #IDEA: Should we add the possibility of having a "color" column specifying 
  #the color of each cytoband? If present, it would take precendence over 
  #"gieStain". Maybe we could also add a separate parameter to specify the 
  #column name of the "gieStain" or color info... but That would grow the 
  #function quite a lot.
  #If the cytobands object does not have the "gieStain" attribute, 
  if(!("gieStain" %in% colnames(mcols(cyto)))) {
    warning("No 'gieStain' column found in cytobands. Using 'gpos50' (gray) for all of them")
    mcols(cyto) <- cbind(mcols(cyto), gieStain="gpos50")
  }
  
  #filter out the cytobands out of our chromosomes
  cyto <- filterChromosomes(cyto, keep.chr = karyoplot$chromosomes)
  
  
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
  ccf <- karyoplot$coord.change.function
  pp <- karyoplot$plot.params
  mids <- karyoplot$ideogram.mid
  
  color.table <- getCytobandColors(color.table, color.schema)
 
  
  
  #And plot them
  ybottom <- mids(as.character(seqnames(cyto)))
  ytop <- mids(as.character(seqnames(cyto)))
    
  xleft <- ccf(x=start(cyto), chr=as.character(seqnames(cyto)))$x
  xright <- ccf(x=end(cyto), chr=as.character(seqnames(cyto)))$x
    
  col <- color.table[as.character(cyto$gieStain)]
  
  if(karyoplot$zoom==TRUE) {
    if(clipping==TRUE) {
      #get the plot coordinates of the cytobands drawing area
      clip.xleft <- ccf(x=start(karyoplot$plot.region), chr=as.character(seqnames(karyoplot$plot.region)))$x
      clip.xright <- ccf(x=end(karyoplot$plot.region), chr=as.character(seqnames(karyoplot$plot.region)))$x
      clip.ybottom <- ybottom - 10 #add a small margin to allow for the width of the lines
      clip.ytop <- ytop + 10
      graphics::clip(x1 = clip.xleft, x2 = clip.xright, y1 = clip.ybottom, y2=clip.ytop)
    }
  }
  graphics::segments(x0 = xleft, x1=xright, y0=ybottom, y1=ytop, col=col, lwd=lwd, lend=lend)      

  
  invisible(karyoplot)
}

