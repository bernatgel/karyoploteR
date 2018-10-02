#' kpAddCytobands
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
## @usage kpAddCytobands(karyoplot, color.table=NULL, add.cytobands.names=FALSE, add.base.numbers=FALSE, ...)
#' @usage kpAddCytobands(karyoplot, color.table=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param color.table  (named character vector) a table specifying the colors to plot the cytobands. If NULL, it gets the colors calling \code{getCytobandColors}. (defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, cytoband representation will be not drawn out of the drawing are (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the cytobands representation may overflow into the margins of the plot. (defaults to TRUE)
## @param add.cytobands.names  (boolean) whether to add or not the cytoband names to the plot. (defaults to FALSE)
## @param add.base.numbers  (boolean) whether to add the base numbers to the plot. (defaults to FALSE)
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
#' kpAddCytobands(kp)
#'  
#' @export kpAddCytobands
#' 



kpAddCytobands <- function(karyoplot, color.table=NULL, clipping=TRUE, ...) {
  
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
    
  color.table <- getCytobandColors(color.table, ...)
  border <- ifelse("border" %in% names(color.table), color.table["border"], "black")
 
  #And plot them
  ybottom <- mids(as.character(seqnames(cyto))) - pp$ideogramheight/2
  ytop <- mids(as.character(seqnames(cyto))) + pp$ideogramheight/2
    
  xleft <- ccf(x=start(cyto), chr=as.character(seqnames(cyto)))$x
  xright <- ccf(x=end(cyto), chr=as.character(seqnames(cyto)))$x
    
  #col <- do.call(c, color.table[cyto$gieStain])
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
  graphics::rect(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop, col=col, border=border)      

  
  invisible(karyoplot)
}

