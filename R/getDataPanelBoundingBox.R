#' getDataPanelBoundingBox
#' 
#' @description 
#' 
#' Return the bounding box of a data panel in plot coordinates
#' 
#' @details 
#' 
#' Given a KaryoPlot object and a data.panel name, return the region 
#' where the data panel is placed. The returned values are in plot
#' coordinates, that is, the coord.change.function has been applied.
#' 
#' @note A user is not expected to need this function. It is mainly
#' used by plotting functions, specially when clipping the plot in 
#' a zoomed region.
#' 
#' @usage getDataPanelBoundingBox(karyoplot, data.panel)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param data.panel   a valid data panel name (i.e. 1, 2, "ideogram", "all")
#' 
#' @return
#' Returns a list with four elements (x0, x1, y0 and y1), each of them an 
#' integer with the plot coordinates for the data.panel
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpDataBackground}}
#' 
#' @examples
#'
#' kp <- plotKaryotype(plot.type=2)
#' dp1 <- getDataPanelBoundingBox(kp, 1)
#'  
#' @export getDataPanelBoundingBox
#' 


getDataPanelBoundingBox <- function(karyoplot, data.panel) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(is.null(karyoplot$plot.params[[paste0("data", data.panel, "min")]])) stop("data panel ", data.panel, " is not available in the current plot.type")
  
  x0 <- karyoplot$coord.change.function(chr=as.character(seqnames(karyoplot$plot.region)), x=start(karyoplot$plot.region), data.panel=data.panel)$x
  x1 <- karyoplot$coord.change.function(chr=as.character(seqnames(karyoplot$plot.region)), x=end(karyoplot$plot.region), data.panel=data.panel)$x
  y0 <- karyoplot$coord.change.function(chr=as.character(seqnames(karyoplot$plot.region)), y=rep(karyoplot$plot.params[[paste0("data", data.panel, "min")]], length(karyoplot$plot.region)), data.panel=data.panel)$y
  y1 <- karyoplot$coord.change.function(chr=as.character(seqnames(karyoplot$plot.region)), y=rep(karyoplot$plot.params[[paste0("data", data.panel, "max")]], length(karyoplot$plot.region)), data.panel=data.panel)$y
  return(list(xleft=xleft, xright=xright, ybottom=ybottom, ytop=ytop))
}