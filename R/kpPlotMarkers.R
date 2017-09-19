#' kpPlotMarkers
#' 
#' @description 
#' 
#' Plots markers on the genome as a line with a label on top.
#' 
#' @details 
#'  
#' This function plots markers on the genome. It implements an interative
#' algorithm to avoid ovelapping between the labels of different markers. Since
#' labels might be plotted in a different position than the original points, 
#' a line with three parts (a vertical, a diagonal and another vertical) is 
#' plotted to link the label with the original position. It is possible to plot
#' labels in horizontal or vertical text and to specify different colors for the
#' marker line and label.
#' 
#' @note The iterative algorithm is not guaranteed to suceed and might end up
#' with overlapping labels if labels are too dense or if too few iterations 
#' allowed. With many markers, the algorithm might be slow.
#'
#' @usage kpPlotMarkers(karyoplot, data=NULL, chr=NULL, x=NULL, y=0.75, labels=NULL, 
#'                      adjust.label.position=TRUE, label.margin=0.001, max.iter=150, label.dist=0.001,
#'                      marker.parts = c(0.8,0.1, 0.1), text.orientation ="vertical",
#'                      ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, 
#'                      line.color="black", label.color="black",
#'                      pos=NULL, srt=NULL, offset=NULL, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the data. If \code{data} is present, \code{chr} will be set to \code{seqnames(data)} and \code{x} to the midpoints of the rages in data. If no parameter \code{y} is specified and \code{data} has a column named \code{y} or \code{value} this column will be used to define the \code{y} value of each data point. (defaults to NULL)
#' @param chr    (a charecter vector) A vector of chromosome names specifying the chromosomes of the data points. If \code{data} is not NULL, \code{chr} is ignored. (defaults to NULL)
#' @param x    (a numeric vector) A numeric vector with the positions (in base pairs) of the data points in the chromosomes. If \code{data} is not NULL, \code{x} is ignored. (defaults to NULL)
#' @param y    (a numeric vector) A numeric vector with the values of the data points. If \code{y} is not NULL, it is used instead of any data column in \code{data}. (defaults to 0.75)
#' @param labels    (a character vector) The labels to be plotted. (defaults to NULL)
#' @param adjust.label.position (logical) whether to adjust the label positions to avoid label overlapping (defaults to TRUE)
#' @param label.margin  (numeric) The vertical margin to leave between the ens of the marker line and the marker label. In plot coordinates. (defaults to 0.001)
#' @param max.iter  (numeric) The maximum number of iterations in the iterative algorithm to adjust the label positioning. (defaults to 150)
#' @param label.dist  (numeric) The minimum distance between labels to consider them as non-overlapping (defaults to 0.001)
#' @param marker.parts   (numeruc vector of three elements) The portion of the distance between 0 and y to be filled with a: vertical, diagonal of vertical part of the marker line. (defaults to c(0.8,0.1,0.1), long vertical stem, small diagonal and small vertical on top)
#' @param text.orientation  ("vertical" or "horizontal) How should the text be plotted. Forced values of srt and pos take precedence. (defaults to "vertical")
#' @param ymin    (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL)
#' @param ymax    (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param line.color  (color) The color of marker line. (defaults to "black")
#' @param label.color  (color) The color of the label (defaults to "black")
#' @param pos (1,2,3,4) The standard pos graphical parameter. If NULL, it's automatically set depending on "text.orientation". (defaults to NULL)
#' @param srt  (numeric) The standard srt graphical parameter. If NULL, it's automatically set depending on "text.orientation". (defaults to NULL)
#' @param offset  (numeric) The standard offset graphical parameter. If NULL, it's automatically set depending on "text.orientation". (defaults to NULL)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#'
#' @return
#' 
#' Returns the original karyoplot object with the data computed (adjusted label positioning) stored at \code{karyoplot$latest.plot}
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{kpLines}}, \code{\link{kpText}}
#' 
#' @examples
#'  
#' 
#' data <- toGRanges(data.frame(c("chr1", "chr1", "chr1"), c(20e6, 21e6, 22e6), c(20.01e6, 21.01e6, 22.01e6), labels=c("GeneA", "GeneB", "GeneC")))
#' 
#' kp <- plotKaryotype("hg19", plot.type=1, chromosomes = "chr1", main="Default markers")
#' kpPlotMarkers(kp, data)
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes = "chr1", main="Markers Horizontal")
#' kpPlotMarkers(kp, data, text.orientation = "horizontal")
#' kpPlotMarkers(kp, data, text.orientation = "horizontal", label.dist = 0.02, data.panel=2)
#' 
#' kp <- plotKaryotype("hg19", plot.type=2, chromosomes = "chr1", main="Different Marker parts")
#' kpPlotMarkers(kp, data, text.orientation = "horizontal", marker.parts=c(0, 1, 0), line.color="red")
#' kpPlotMarkers(kp, data, text.orientation = "horizontal", marker.parts=c(0.1, 0.2, 0.4), label.dist = 0.02, data.panel=2, label.color="blue")
#' 
#' 
#'  
#' @export kpPlotMarkers
#' @importFrom digest digest
#' @importFrom graphics strwidth

#TODO: Adjust the positioning algorithm when zoom is active so labels do not fall out of the plot.region

kpPlotMarkers <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=0.75, labels=NULL, 
                   adjust.label.position=TRUE, label.margin=0.001, max.iter=150, label.dist=0.001,
                   marker.parts = c(0.8,0.1, 0.1), text.orientation ="vertical",
                   ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, 
                   line.color="black", label.color="black",
                   pos=NULL, srt=NULL, offset=NULL, clipping=TRUE, ...) {
  
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  text.orientation <- match.arg(text.orientation, c("horizontal", "vertical"))
  
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  if(!is.null(data)) {
    if("labels" %in% names(mcols(data))) {
      if(is.null(labels) & (!is.null(data) & (is(mcols(data)[,"labels"], "factor") | is(mcols(data)[,"labels"], "character")))) {
        labels <- as.character(mcols(data)[,"labels"])
      }
      if(is.null(labels) & (!is.null(data) & (is(mcols(data)[,1], "factor") | is(mcols(data)[,1], "character")))) {
        labels <- as.character(mcols(data)[,1])
      }
    }
  }
  if(is.null(labels)) {
    stop("labels is NULL and no valid labels found in data")
  }
  
  pp <- prepareParameters2("kpPlotMarkers", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y,
                           ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  #Filter the data so we only have the data in the plot.region
  valid <- which(overlapsAny(GRanges(pp$chr, IRanges(pp$x, pp$x)), karyoplot$plot.region))
  
  if(length(valid)==0) { #If there are no markers in the plot region, just return
    karyoplot$latest.plot <- list(funct="kpPlotMarkers", computed.values=list(label.position=integer(0)))
    invisible(karyoplot)
  }
  
  valid.chr <- pp$chr[valid]
  valid.x <- pp$x[valid]
  valid.y <- pp$y[valid]
  valid.labels <- labels[valid]
  
  #Sort everything by chr and start
  ord <- order(valid.chr, valid.x)
  valid.chr <- valid.chr[ord]
  valid.x <- valid.x[ord]
  valid.y <- valid.y[ord]
  valid.labels <- valid.labels[ord]
  
  #Transform to plot coordinates  
  xplot <- ccf(chr=valid.chr, x=valid.x, data.panel=data.panel)$x
  yplot <- ccf(chr=valid.chr, y=valid.y, data.panel=data.panel)$y
  
  if(adjust.label.position==TRUE) {
      
    #Now, move the text as needed (only horizontally) to avoid overlapping
    xlabels <- xplot
    
    for(current.chr in karyoplot$chromosomes) {
      in.chr <- valid.chr==current.chr
      
      if(any(in.chr)) {
        
        xp <- xlabels[in.chr]
        chr.labels <- valid.labels[in.chr]
        
        if(text.orientation == "horizontal") {
          sw <- strwidth(chr.labels, units = "user") + label.dist
        } else {
          sw <- strwidth(rep("M", length(chr.labels)), units = "user") + label.dist
        }
        delta <- strwidth("a", units = "user")/4
        
        #to be sure that we do not enter a loop, store a hash of the all old positions
        old.pos <- digest(xp)
    
        for(iter in seq_len(max.iter)) {
          for(i in seq_len(length(chr.labels)-1)) {
            #Test for overlap with the next label
            if((xp[i]+sw[i]/2) > (xp[i+1]-sw[i+1]/2)) {
              #Decide what to move
              
              #If it's the first label
              if(i==1) {
                #if it will not fall out of the chromosome, move the label to the left
                if((xp[i]-sw[i]/2-delta) > ccf(chr=current.chr, x=0, data.panel=data.panel)$x) {
                  xp[i] <- xp[i]-delta
                } else {
                  #if it will not fall out of the chromosome by the right side
                  if(xp[i+1]+delta < ccf(chr=current.chr, x=karyoplot$chromosome.lengths[current.chr], data.panel=data.panel)$x) {
                    #move the next label to right
                    xp[i+1] <- xp[i+1]+delta  
                  } else {
                    #do nothing
                  }
                }
              } else { #if it's not the first label
                #if it will not overlap with the previous one, move to the left
                if(!((xp[i]-sw[i]/2-delta) < (xp[i-1]+sw[i-1]/2))) { 
                  xp[i] <- xp[i]-delta
                } else {
                  #if it will not fall out of the chromosome by the right side
                  if(xp[i+1]+delta < ccf(chr=current.chr, x=karyoplot$chromosome.lengths[current.chr], data.panel=data.panel)$x) {
                    #move the next label to right
                    xp[i+1] <- xp[i+1]+delta  
                  } else {
                    #do nothing
                  }
                }
              }
            }
          }
          dd <- digest(xp)
          if(dd %in% old.pos) { #We have either found a valid position or entered a loop
            break; 
          } else {
            old.pos <- c(old.pos, dd)
          }
        }
        
        xlabels[in.chr] <- xp
      }
    }
  } else { #if we do not reposition the labels, the label position is the same as the original one
    xlabels <- xplot
  }
  #Finished repositioning labels

  if(karyoplot$zoom==TRUE) {
    if(clipping==TRUE) {
      dpbb <- karyoplot$getDataPanelBoundingBox(data.panel)
      graphics::clip(x1 = dpbb$xleft, x2 = dpbb$xright, y1 = dpbb$ybottom, y2=dpbb$ytop)
    }
  }
  
  #Plot the labels
  #detect if the data.panel is an inverted one
  is.inverted <- ccf(chr=karyoplot$chromosomes[1], y=0, data.panel=data.panel)$y > ccf(chr=karyoplot$chromosomes[1], y=1, data.panel=data.panel)$y
  if(text.orientation == "horizontal") {
    if(is.inverted==FALSE) {
      pos <- ifelse(is.null(pos), 3, pos)
      label.y <- yplot+label.margin
    } else {
      pos <- ifelse(is.null(pos), 1, pos)
      label.y <- yplot-label.margin
    }
    graphics::text(x=xlabels, y=label.y, labels=valid.labels, pos=pos, col=label.color, ...)
  } else {
    srt <- ifelse(is.null(srt), 90, srt)
    offset <- ifelse(is.null(offset), 0, offset)
    if(is.inverted==FALSE) {
      pos <- ifelse(is.null(pos), 4, pos)
      label.y <- yplot+label.margin
    } else {
      pos <- ifelse(is.null(pos), 2, pos)
      label.y <- yplot-label.margin
    }
    graphics::text(x=xlabels, y=label.y, labels=valid.labels, pos=pos, srt=srt, offset=offset, col=label.color, ...)

  }
  #and plot the markers
  if(marker.parts[2]==0 & adjust.label.position==TRUE) { #If no diagonal part is allowed, warn the user
    if(!isTRUE(all.equal(valid.x, xlabels))) {
      warning("Markers labels have been repositioned, but the diagonal part of the marker line is set to 0. For better results set the second value of 'marker.parts' to something different than 0.")
    }
  }
  
  #divide the vertical space according to the proportions of the parts
  my <- marker.parts*1/sum(marker.parts)
  
  #and plot them
    #Determine the y of the value y=0 to start the marker there
      pp2 <- prepareParameters2("kpPlotMarkers", karyoplot=karyoplot, data=NULL, chr=pp$chr, x=pp$x, y=0, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
      y0 <-  ccf(chr=pp2$chr[valid][ord], y=pp2$y[valid][ord], data.panel=data.panel)$y
    #markers
    if(marker.parts[1]>0) {
      graphics::segments(x0 = xplot, x1=xplot, y0 = y0, y1=y0+my[1]*(yplot-y0), col=line.color, ... )
    }
    if(marker.parts[2]>0) {
      graphics::segments(x0 = xplot, x1=xlabels, y0 = y0+my[1]*(yplot-y0), y1=y0+(my[1]+my[2])*(yplot-y0), col=line.color, ...)
    }
    if(marker.parts[3]>0) {
      graphics::segments(x0 = xlabels, x1=xlabels, y0 = yplot - my[3]*(yplot-y0), y1=yplot, col=line.color, ...)
    }
    
  karyoplot$latest.plot <- list(funct="kpPlotMarkers", computed.values=list(label.position=xlabels))
    
  invisible(karyoplot)
}
