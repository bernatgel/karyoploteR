#' kpChars
#'
#' @description
#'
#' Plot bars along the genome
#'
#' @details
#'
#' \code{kpChars} plots bars (rectangles) along the genome. It is very similar to
#' \code{\link{kpRect}} except that if \code{y0} is missing, it's automatically set
#' to \code{ymin} so all bars start from the base of the plotting region.
#'
#' @usage kpChars(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y1=NULL, y0=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, clipping=TRUE, ...)
#'
#' @inheritParams kpRect
#'
#' @return
#'
#' Returns the original karyoplot object, unchanged.
#'
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpRect}}, \code{\link{kpLines}}
#'
#' @examples
#'

#'
#' @export kpChars

#Plot single characters from a nucleic acid sequence expanding the full rectangle defined by x0-x1 and y0-y1
#Can plot multiple chars if any of the arguments is a vector and everything else will be recycled. Right? Quite low level. Allows for logos and SNPs

kpChars <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL,
                   chars=NULL, plot.letters=TRUE, plot.rectangles=FALSE,
                   letter.colors = "IGV", letter.borders = NA,
                   rectangle.colors="IGV", rectangle.borders=NA,
                   ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL,
                   clipping=TRUE, ...) {

  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpChars: 'karyoplot' must be a valid 'KaryoPlot' object"))

  #TODO: Validation
  if(length(chars)==0 || any(is.na(chars)) || any(is.null(chars))) stop("In kpChars: you must provide at least one character in 'chars'")
  if(!all(is.character(chars))) stop("In kpChars: 'chars' must be a character or character vector")
  #if chars is a single character string with more than one character, split it into individual chars
  if(length(chars)==1 & nchar(chars)>1) {
    chars <- strsplit(chars, "")[[1]]
  }
  if(!all(unlist(lapply(chars, length))==1)) stop("In kpChars: all elements in 'chars' must be single characters. To plot multiple characters you must provide a a single character string or a vector of individual characters (i.e. c('A', 'T', 'T', 'G'))")

  #If x0 is not specified with any of the valid methods, error
  if(is.null(x0)) {
    if(is.null(data)) {
      stop("In kpChars: the start position for each character is needed. You must specify it using 'x0' or 'data'.")
    } else { #If data exist, then x0 will be the start of each range
      x0=GenomicRanges::start(data)
    }
  }

  #If x1 is not specified with any of the valid methods, set it to the x0+1, so every char occupies exactly one base
  if(is.null(x1)) {
    if(is.null(data)) {
      x1=x0+1
    } else { #If data exist, the x1 will be the end of each range plus one (so the letters end just before the start of the end base)
      x1=GenomicRanges::end(data)+1
    }
  }

  #If y0 is not specified with any of the valid methods, set it to the min of the data.panel
  if(is.null(y0)) {
    if(is.null(data)) {
      y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
    } else {
      if(!("y0" %in% names(mcols(data)))) {
        y0=karyoplot$plot.params[[paste0("data", data.panel, "min")]]
      }
    }
  }

  #If y1 is not specified with any of the valid methods, set it to the max of the data.panel
  if(is.null(y1)) {
    if(is.null(data)) {
      y1=karyoplot$plot.params[[paste0("data", data.panel, "max")]]
    } else {
      if(!("y1" %in% names(mcols(data)))) {
        y1=karyoplot$plot.params[[paste0("data", data.panel, "max")]]
      }
    }
  }

  #Now check that x0 x1 y0 and y1 do not containg any NA and are numbers


  pp <- prepareParameters4("kpChar", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1,
                           y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1,
                           data.panel=data.panel, ...)
 message("X0: ", paste(x0, collapse = ","))
 message("X1: ", paste(x1, collapse = ","))
 message("Y0: ", paste(y0, collapse = ","))
 message("Y1: ", paste(y1, collapse = ","))

  #Get the color schemas
  #and get the color for each char and rectangle
  if(plot.rectangles==TRUE) {
    rect.col.schema <-  getColorSchemas()$sequences$schemas[[rectangle.colors]]

    #TODO
    #if(!is.null(rectangle.colors)
    rect.cols <- rect.col.schema[chars]
    rect.brd <- NA
    #TODO: border colors
  }
  if(plot.letters==TRUE) {
    #If letter.colors is an array of colors, use them to plot the letters
    if(all(is.color(letter.colors) | is.na(letter.colors))) {
      letter.cols <- letter.colors
      #TODO: Check if it's only a color vector or a named color vector we can use as a schema
    } else {
      #Is it the name of a schema?
      if(!(letter.colors %in% names(getColorSchemas()$sequences$schemas))) {
        warning("No valid colors or color schema found for letters. Using black as default. You can provide it with 'letter.colors'.")
        letter.cols <- "black"
      } else {
        letter.col.schema <- getColorSchemas()$sequences$schemas[[letter.colors]]
        letter.cols <- letter.col.schema[chars]
        #If the color schema has a default color, change all NAs to the default color
        if("default" %in% names(letter.col.schema)) {
          letter.cols[is.na(letter.cols)] <- letter.col.schema["default"]
        }
      }
      message(letter.cols)
    }
    #TODO: repeat the same for borders

    letter.brd <- NA
    #TODO: Use defaults for those with NULL or NA
  }
  #Repeat the same for rectangle colors and borders


  #Plot the rectangles
  if(plot.rectangles==TRUE) {
    kpRect(karyoplot, chr=pp$chr, x0=pp$x0, x1=pp$x1, y0=pp$y0, y1=pp$y1,
           col=rect.cols, border=rect.brd, data.panel = data.panel,
           clipping = clipping, ...)
  }

  if(plot.letters) {

    #Call kpPolygon once per character to plot
    for(nc in seq_along(chars)) {
      cpol <- char.polygons[[chars[nc]]]
      #Since we have already used prepareParameters, we set r's and ymin and ymax to 0 and 1.
      kpPolygon(karyoplot, chr=pp$chr[nc], x=pp$x0[nc]+cpol$x*(pp$x1[nc]-pp$x0[nc]), y=pp$y0[nc]+cpol$y*(pp$y1[nc]-pp$y0[nc]),
                col=letter.cols[nc], border=letter.brd,
                r0=0, r1=1, ymin = 0, ymax=1, data.panel = data.panel, clipping=clipping )
      #TODO: set correct color for the border
    }
  }
  #TODO: store the info in latest plot
  #including: letter and rectangle colors
  invisible(karyoplot)
}






#
#
#
#
# #Plot a sequence. Uses kpChar for the actual plot
# #No
# kpPlotSequence <- function(karyoplot, data=NULL, chr=NULL, x=NULL, sequence=NULL, height=1.2,
#                            plot.letters=TRUE, plot.rectangles=FALSE,
#                            letter.colors = "IGV", letter.borders = NA,
#                            rectangle.colors="IGV", rectangle.borders=NA,
#                            data.panel = 1, r0 = NULL, r1 = NULL, clipping = TRUE, ...) {
#
#   y.mult <- 1
#
#   #convert the sequence to an array of characters
#   seq.chars <- strsplit(sequence, "")[[1]]
#   for(nc in seq_along(seq.chars)) {
#     char <- seq.chars[nc]
#
#   }
#
#
# }
