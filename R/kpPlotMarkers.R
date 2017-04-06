#' kpPlotMarkers
#' 
#' @description 
#' 
#' Plots the text given in \code{labels} at the positions defined by chr, x and y along the genome.
#' 
#' @details 
#'  
#' This is one of the functions from karyoploteR implementing the adaptation to the genome context 
#' of basic plot functions
#' from R base graphics. Given a set of positions on the genome (chromosome and base), a 
#' value (y) for each of them and a label, it plots the label at the position specified by the 
#' data point. Data can be provided via a \code{GRanges} object (\code{data}), independent
#' parameters for chr, x and y or a combination of both. A number of parameters can be used 
#' to define exactly where 
#' and how the text is drawn. In addition, via the ellipsis operator (\code{...}), \code{kpPlotMarkers}
#' accepts any parameter valid for \code{text} (e.g. \code{cex}, \code{col}, ...)
#'
#' @usage kpPlotMarkers(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, labels=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...)
#' 
#' @param labels    (a character vector) The labels to be plotted. (defaults to NULL)
#' @inheritParams kpPoints
#'
#' @return
#' 
#' Returns the original karyoplot object, unchanged. 
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpPoints}}
#' @seealso \code{\link{kpPlotRegions}}
#' 
#' @examples
#'  
#' set.seed(1000)
#' data.points <- sort(createRandomRegions(nregions=500, mask=NA))
#' mcols(data.points) <- data.frame(y=runif(500, min=0, max=1))
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes=c("chr1", "chr2"))
#'   kpDataBackground(kp, data.panel=1)
#'   kpDataBackground(kp, data.panel=2)
#' 
#'   kpLines(kp, data=data.points, col="red")
#' 
#'   #Three ways of specifying the exact same data.points
#'   kpPoints(kp, data=data.points)
#'   kpPoints(kp, data=data.points, y=data.points$y, pch=16, col="#CCCCFF", cex=0.6)
#'   kpPoints(kp, chr=as.character(seqnames(data.points)), 
#'            x=(start(data.points)+end(data.points))/2,
#'            y=data.points$y, pch=".", col="black", cex=1)
#' 
#'   #plotting in the data.panel=2 and using r0 and r1, ymin and ymax
#'   kpLines(kp, data=data.points, col="red", r0=0, r1=0.3, data.panel=2)
#'   kpPlotMarkers(kp, data=data.points, labels=as.character(1:500), r0=0, r1=0.3, data.panel=2, pch=".", cex=3)
#' 
#'   kpLines(kp, data=data.points, col="blue", r0=0.4, r1=0.7, data.panel=2)
#'   kpLines(kp, data=data.points, col="blue", y=-1*(data.points$y), ymin=-1, ymax=0, r0=0.7, r1=1, data.panel=2)
#'   #It is also possible to "flip" the data by giving an r0 > r1
#'   kpPoints(kp, data=data.points, col="red", y=(data.points$y), r0=1, r1=0.7, data.panel=2, pch=".", cex=2)  
#' 
#' 
#'  
#' @export kpPlotMarkers
#' @importFrom digest digest


kpPlotMarkers <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=0.75, labels=NULL, 
                   adjust.label.position=TRUE, label.margin=0.001, max.iter=150,
                   marker.parts = c(0,0.8, 0.2), text.orientation ="horizontal",
                   ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  text.orientation <- match.arg(text.orientation, c("horizontal", "vertical"))
  
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  if(is.null(labels) & (!is.null(data) & (is(mcols(data)[,"labels"], "factor") | is(mcols(data)[,"labels"], "character")))) {
    labels <- as.character(mcols(data)[,1])
  }
  if(is.null(labels) & (!is.null(data) & (is(mcols(data)[,1], "factor") | is(mcols(data)[,1], "character")))) {
    labels <- as.character(mcols(data)[,1])
  }
  if(is.null(labels)) {
    stop("labels is NULL and no valid labels found in data")
  }
  
  pp <- prepareParameters2("kpPlotMarkers", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y,
                           ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
    
  xplot <- ccf(chr=pp$chr, x=pp$x, data.panel=data.panel)$x
  yplot <- ccf(chr=pp$chr, y=pp$y, data.panel=data.panel)$y
  
  if(adjust.label.position==TRUE) {
      
    #Now, move the text as needed (only horizontally) to avoid overlapping
    xp <- xplot
    if(text.orientation == "horizontal") {
      sw <- strwidth(labels, units = "user") + 0.001
    } else {
      sw <- strwidth(rep("aa", length(labels)), units = "user") + 0.001
    }
    delta <- strwidth("a", units = "user")/4
    
    #to be sure that we do not enter a loop, store a hash of the all old positions
    old.pos <- digest(xp)

    for(iter in seq_len(max.iter)) {
      for(i in seq_len(length(labels)-1)) {
        #Test for overlap with the next label
        if((xp[i]+sw[i]/2) > (xp[i+1]-sw[i+1]/2)) {
          #Decide what to move
          
          #If it's the first label
          if(i==1) {
            #if it will not fall out of the chromosome, move the label to the left
            if((xp[i]-sw[i]/2-delta) > ccf(chr=pp$chr[i], x=0, data.panel=data.panel)$x) {
              message("(",i,") ", i, " left 1")
              xp[i] <- xp[i]-delta
            } else { #else, move the next label to to right
              message("(",i,") ", i+1, " right 1")
              xp[i+1] <- xp[i+1]+delta
            }
          } else { #if it's not the first label
            #if it will not overlap with the previous one, move to the left
            if(!((xp[i]-sw[i]/2-delta) < (xp[i-1]+sw[i-1]/2))) { 
              message("(",i,") ", i, " left 2")
              xp[i] <- xp[i]-delta
            } else {
              #move the next label to right
              message("(",i,") ", i+1, " right 2")
              xp[i+1] <- xp[i+1]+delta  
            }
          }
        }
      }
      print(xp)
      dd <- digest(xp)
      if(dd %in% old.pos) { #We have either found a valid position or entered a loop
        break; 
      } else {
        old.pos <- c(old.pos, dd)
      }
    }
    xlabels <- xp
  } else { #if we do not reposition the labels, the label position is the same as the original one
    xlabels <- xplot
  }
  #Finished repositioning labels

  
  #Plot the labels
  if(text.orientation == "horizontal") {
    graphics::text(x=xlabels, y=yplot+label.margin, labels=labels, pos=3, ...)
  } else {
    graphics::text(x=xlabels, y=yplot+label.margin, labels=labels, pos=4, srt=90, offset=0, ...)
  }
  #and plot the markers
  m1 <- marker.parts[1]
  m2 <- marker.parts[2]
  m3 <- marker.parts[3]
  if(m1==0 & m3==0 & adjust.label.position==TRUE) { #If no diagonal part is allowed, make the whole marker diagonal
    m1 <- m2 
    m2 <- 0
  }
  
  #Distribute the x distance between the first and third part of the marker
  if(m1==0 & m3!=0) {m2x <- xplot}
  if(m1!=0 & m3!=0) {m2x <- xplot + (xlabels - xplot)/2}
  if(m1!=0 & m3==0) {m2x <- xlabels}
  
  #divide the vertacal space according to the proportions of the parts
  my <- marker.parts*1/sum(marker.parts)
  
  #and plot them
    #Determine the y of the value y=0 to start the marker there
    pp2 <- prepareParameters2("kpPlotMarkers", karyoplot=karyoplot, data=data, chr=chr, x=x, y=0, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
    y0 <-  ccf(chr=pp2$chr, y=pp2$y, data.panel=data.panel)$y
    #marker 1
    if(m1>0) {
      graphics::segments(x0 = xplot, x1=m2x, y0 = y0, y1=y0+my[1]*(yplot-y0), ... )
    }
    if(m2>0) {
      graphics::segments(x0 = m2x, x1=m2x, y0 = y0+my[1]*(yplot-y0), y1=y0+(my[1]+my[2])*(yplot-y0), ...)
    }
    if(m3>0) {
      graphics::segments(x0 = m2x, x1=xlabels, y0 = yplot - my[3]*(yplot-y0), y1=yplot, ...)
    }
    
  invisible(karyoplot)
}
