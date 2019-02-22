#' getVariantsColors
#' 
#' @description 
#' 
#' Given the reference and alternative for a set of variants, assigns a color 
#' to each of them
#' 
#' @details 
#'  
#' The function creates an nucleotide substitution identifier with for each 
#' variant and uses it to query the color.table lookup table. If color.table
#' is NULL, a color.table based in the selected color.schema is used. All 
#' unkwonwn nucleotide substitutions are assigned a gray color. Color table 
#' needs to have entries for C>A, C>G, C>T, T>A, T>C and T>G (and optionally 
#' "others"), since other changes can be reverse complemented to these.
#' 
#' @usage getVariantsColors(ref, alt, color.table=NULL, color.schema=c("cell21breast"))
#' 
#' @param ref (character vector) The reference nucleotides of the variants. It has to have the same length as \code{alt}.
#' @param alt (character vector) The alternative nucleotides of the variants. It has to have the same length as \code{ref}
#' @param color.table   (named character vector) if present, its used to assign colors to the nucleotide substitutions.
#' @param color.schema  (character) The name of the color schema to use: \code{cell21breast} (the color schema used in "Mutational Processes Molding the Genomes of 21 Breast Cancers" by S. Nik-Zainal, Cell, 2012). (defaults to \code{cell21breast})
#' 
#' @return
#' a named character vector with the colors associated to each variant
#'   
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpPlotRainfall}}
#' 
#' @examples
#'  
#' ref <- c("A", "A", "C", "T", "G", "A")
#' alt <- c("G", "C", "T", "A", "A", "-")
#' getVariantsColors(ref, alt)
#' 
#' col.table <- c("C>A"="#FF0000", "C>G"="#000000", "C>T"="#00FF00", "T>A"="#0000FF", "T>C"="#BB00BB", "T>G"="#00BBBB", "other"="#888888")
#' getVariantsColors(ref, alt, col.table)
#' 
#' 
#' @export getVariantsColors
#' 


getVariantsColors <- function(ref, alt, color.table=NULL, color.schema=c("cell21breast")) {
  
  if(!methods::is(ref, "character")) stop(paste0("In getVariantsColors: 'ref' must be a valid character object"))
  if(!methods::is(alt, "character")) stop(paste0("In getVariantsColors: 'alt' must be a valid character object"))
  if(length(ref) != length(alt)) stop(paste0("In getVariantsColors: 'ref' and 'alt' must have the same length"))
  
  if(is.null(color.table)) {
    color.schema <- match.arg(color.schema)
    
    if(color.schema %in% names(.karyoploter.colors$variants$schemas)) {
      color.table <- .karyoploter.colors$variants$schemas[[color.schema]]
    } else {
      stop("Unknown color.schema. Available schemas for variants are: ", paste0(names(.karyoploter.colors$variants$schemas), collapse = ", "))
    }
  }
    
  comp <- c(G="C", A="T", C="G", T="A")
  
  #complement if necessary
  to.comp <- which(ref=="A" | ref=="G")
  ref[to.comp] <- comp[ref[to.comp]]
  alt[to.comp] <- comp[alt[to.comp]]
  
  var.cols <- color.table[paste0(ref, ">", alt)]
  if(!is.null(color.table["other"])) {
    var.cols[which(is.na(var.cols))] <- color.table["other"]
  } else {
    var.cols[which(is.na(var.cols))] <- "#888888"
  }
  
  return(var.cols)
}


