#' getCytobandColors
#' 
#' @description 
#' 
#' Returns a named character vector with the colors of associated with the cytoband names
#' 
#' @details 
#'  
#' The function returns a named character vector with the colors of associated with the 
#' cytoband names. Two color schemas are available: circos (which copies the colors used by 
#' Circos) and biovizbase (that gets the cytoband colors from the biovizBase Bioconductor 
#' package). If a color.table is given, it is returned untouched.
#' 
#' @usage getCytobandColors(color.table=NULL, color.schema=c("circos", "biovizbase", "only.centromeres")) 
#' 
#' @param color.table   (named character vector) if present, it's returned as-is. Useful to specify your own color.tables.
#' @param color.schema  (character) The name of the color schema to use: \code{circos}, \code{biovizBase}, \code{only.centromeres} (everything in gray, except for centromeres in red). (defaults to \code{circos})
#' 
#' @return
#' a named character vector with the colors associated to each cytoband name
#' 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{plotCytobands}}
#' 
#' @examples
#'  
#' getCytobandColors()
#' getCytobandColors(color.schema="biovizbase")
#'  
#' @export getCytobandColors
#' 


getCytobandColors <- function(color.table=NULL, color.schema=c("circos", "biovizbase", "only.centromeres")) {
  
  color.schema <- match.arg(color.schema)
  
  if(is.null(color.table)) {
    if(color.schema=="biovizbase") {
      color.table <- biovizBase::getBioColor("CYTOBAND")
    } else {
      if(color.schema=="circos") {
         color.table <- c(gneg="#FFFFFF",
                         gpos25="#C8C8C8",
                         gpos33="#D2D2D2",
                         gpos50="#C8C8C8",
                         gpos66="#A0A0A0",
                         gpos75="#828282",
                         gpos100="#000000",
                         gpos="#000000",
                         stalk="#647FA4", #repetitive areas
                         acen="#D92F27", #centromeres
                         gvar="#DCDCDC")
      } else {
        if(color.schema=="only.centromeres") {
          color.table <- c(gneg="#C8C8C8",
                           gpos25="#C8C8C8",
                           gpos33="#C8C8C8",
                           gpos50="#C8C8C8",
                           gpos66="#C8C8C8",
                           gpos75="#C8C8C8",
                           gpos100="#C8C8C8",
                           gpos="#C8C8C8",
                           stalk="#C8C8C8", #repetitive areas
                           acen="#D92F27", #centromeres
                           gvar="#C8C8C8")
        }
      }
    }
  }
    
  return(color.table)
}


