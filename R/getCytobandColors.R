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
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpAddCytobands}}
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
  
  if(!is.null(color.table)) {
    #TODO: Should we check if it makes sense or leave it to the user?
    return(color.table)
  } else {
    if(color.schema %in% names(.karyoploter.colors$cytobands$schemas)) {
      return(.karyoploter.colors$cytobands$schemas[[color.schema]])
    } else {
      stop("Unknown color.schema. Available schemas for cytobands are: ", paste0(names(.karyoploter.colors$cytobands$schemas), collapse = ", "))
    }
  }
  #We should never reach this point. 
  stop("Error in getCytobandColors")
}
