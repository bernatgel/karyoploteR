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
#' set.seed(1000)
#' 
#' data <- toGRanges(data.frame(chr="chr1", start=10000000*(0:23), end=10000000*(1:24)))
#' y1 <- ((sin(start(data)) + rnorm(n=24, mean=0, sd=0.1))/5)+0.5
#' y0 <- y1 - rnorm(n=24, mean = 0, sd = 0.15)
#'  
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#' 
#' #We can specify all data values separately. If missing y0, it defaults to ymin
#' kpChars(kp, chr=as.character(seqnames(data)), x0=start(data), x1=end(data), y1=y1,
#'        col="#FFBBBB", border="#EEAAAA")
#' kpLines(kp, data=data, y=y1, col="red")
#' 
#' #or we can provide all data into a single GRanges object
#' mcols(data) <- data.frame(y0=y0, y1=y1)
#' kpChars(kp, data[data$y0>data$y1], col="orange", border="orange", data.panel=2)
#' kpChars(kp, data[data$y0<=data$y1], col="purple", border="purple", data.panel=2)
#' 
#' kpLines(kp, data, y=data$y1, data.panel=2, col="red")
#' kpLines(kp, data, y=data$y0, data.panel=2, col="blue")
#' 
#' kpAxis(kp, data.panel = 1, cex=0.8, numticks = 5, col="#777777")
#' kpAxis(kp, data.panel = 2, cex=0.8, numticks = 5, col="#777777")
#' 
#' @export kpChars

#Plot single characters from a nucleic acid sequence expanding the full rectangle defined by x0-x1 and y0-y1 
#Can plot multiple chars if any of the arguments is a vector and everything else will be recycled. Right? Quite low level. Allows for logos and SNPs

kpChars <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y0=0, y1=1,
                   chars=NULL, plot.letters=TRUE, plot.rectangles=FALSE,
                   letter.colors = "IGV", letter.borders = NA,
                   rectangle.colors="IGV", rectangle.borders=NA,
                   ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, 
                   clipping=TRUE, ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop(paste0("In kpChars: 'karyoplot' must be a valid 'KaryoPlot' object"))

  #TODO: Validation
  if(!length(chars)==0 || is.na(char)) stop("In kpChars: you must provide at least one character in 'chars'")
  if(!all(is.character(chars))) stop("In kpChars: 'chars' must be a character or character vector")
  if(!all(unlist(lapply(chars, length))==1)) stop("In kpChars: all elements in 'chars' must be single characters. To plot multiple characters you must provide a vector of characters (i.e. c('A', 'T', 'T', 'G'))")

  
  pp <- prepareParameters4("kpChar", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1,
                           y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, 
                           data.panel=data.panel, ...)

    
  
  #Get the color schemas
  #TODO: get them only if needed. Not if NA. respect color provided by user, if any
  letter.col.schema <- getColorSchemas()$sequences$schemas[[letter.colors]]
  rect.col.schema <-  getColorSchemas()$sequences$schemas[[rectangle.colors]]
  
  #and get the color for each char
  if(plot.rectangles==TRUE) {
    #TODO
    #if(!is.null(rectangle.colors) 
    rect.cols <- rect.col.schema[char]
    rect.border <- NA
    #TODO: border colors
  }
  if(plot.letters==TRUE) {
    #TODO
    letter.cols <- letter.col.schema[char] #Use defaults for those wit NULL or NA
  }
    
  pp
  #Plot the rectangles
  if(plot.rectangles==TRUE) {
    kpRect(karyoplot, chr=pp$chr, x0=pp$x0, x1=pp$x1, y0=pp$y0, y1=pp$y1, col=rect.cols, border=rect.border, data.panel = data.panel, clipping = clipping), ...)
  }
    
  if(plot.letters) {
    
    #Call kpPolygon once per character to plot
    for(nc in seq_along(chars)) {
      cpol <- char.polygons[[chars[nc]]]
      #Since we have already used prepareParameters, we set r's and ymin and ymax to 0 and 1.
      kpPolygon(karyoplot, chr=pp$chr[nc], x=pp$x0[nc]+cpol$x*(pp$x1[nc]-pp$x0[nc]), y=pp$y0[nc]+cpol$y*(pp$y1[nc]-pp$y0[nc]), 
                col=letter.cols[nc], border="black",
                r0=0, r1=1, ymin = 0, ymax=1, data.panel = data.panel, clipping=clipping )
      #TODO: set correct color for the border
    }
    
  }
  

  
  
  
  
  char.poly <- alphabet[alphabet$group==char,]
  
  
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
  
  
  
  #TODO: we are assuming only one char
  col <- letter.colors
  char.poly <- alphabet[alphabet$group==char,]
  
  

  
  
  invisible(kpRect(karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1,
                   ymin=ymin, ymax=ymax, r0=r0, r1=r1,
                   data.panel=data.panel, clipping=clipping, ...))
  
}





#Plot a sequence. Uses kpChar for the actual plot
#No 
kpPlotSequence <- function(karyoplot, data=NULL, chr=NULL, x=NULL, sequence=NULL, height=1.2, 
                           plot.letters=TRUE, plot.rectangles=FALSE,
                           letter.colors = "IGV", letter.borders = NA,
                           rectangle.colors="IGV", rectangle.borders=NA,
                           data.panel = 1, r0 = NULL, r1 = NULL, clipping = TRUE, ...) {
  
  y.mult <- 1
  
  #convert the sequence to an array of characters
  seq.chars <- strsplit(sequence, "")[[1]]
  for(nc in seq_along(seq.chars)) {
    char <- seq.chars[nc]
    
  }
  
  
}
