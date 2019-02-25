############  Text related functions for karyoploteR   ###############

#' getTextSize
#' 
#' @description 
#' Returns the size of character strings in bases and r's
#' 
#' @details 
#' Small utility function to get the size of text labels in usable units for
#' karyoploteR: bases for the width and r's for the height. The r units are the
#' ones passed to r0 and r1 and take into account that a data panel has always
#' total height of 1 r (from r0=0 to r1=1)
#' 
#' @usage getTextSize(karyoplot, labels, cex=1, data.panel="1")
#' @param karyoplot (KaryoPlot) A KaryoPlot object representing the current plot
#' @param labels  (character) The character strings to measure
#' @param cex (numeric) The cex value used to plot the text (defaults to 1)
#' @param data.panel (data panel identifier) The name of the data panel on which text was plotted. (defaults to "1")
#'
#' @return
#' Returns a list with two elements: width and height. Each of them is a numeric
#' vector of the same length as "labels" with the width in bases of each label
#' and the height in r units of each label.
#' 
#' @examples
#' 
#' pp <- getDefaultPlotParams(plot.type=2)
#' pp$data2height <- 50
#' 
#' kp <- plotKaryotype(chromosomes="chr1", plot.type=2, plot.params=pp)
#'  
#' label <- "Looooooong label"
#' kpText(kp, chr="chr1", x=70e6, y=0.5, labels=label)
#' text.size <- getTextSize(kp, labels=label)
#' kpRect(kp, chr="chr1", x0=70e6-text.size$width/2, x1=70e6+text.size$width/2,
#'            y0=0.5-text.size$height/2, y1=0.5+text.size$height/2)
#'  
#' label <- "SHORT"
#' text.size <- getTextSize(kp, labels=label, cex=3)
#' kpRect(kp, chr="chr1", x0=170e6-text.size$width/2, x1=170e6+text.size$width/2,
#'            y0=0.2-text.size$height/2, y1=0.2+text.size$height/2, col="gold")
#' kpText(kp, chr="chr1", x=170e6, y=0.2, labels=label, cex=3)
#'     
#' label <- c("two_labels", "in a small data.panel=2")
#' kpText(kp, chr="chr1", x=c(100e6, 170e6), y=c(0.4, 0.2), labels=label, cex=0.6, data.panel=2)
#' text.size <- getTextSize(kp, labels=label, cex=0.6, data.panel=2)
#' kpRect(kp, chr="chr1", x0=c(100e6, 170e6)-text.size$width/2, x1=c(100e6, 170e6)+text.size$width/2,
#'            y0=c(0.4, 0.2)-text.size$height/2, y1=c(0.4, 0.2)+text.size$height/2, data.panel=2)
#' 
#'     
#' @export getTextSize
#'


getTextSize <- function(karyoplot, labels, cex=1, data.panel="1") {
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  if(length(labels)==0) return(NULL)
  if(!all(is.character(labels))) stop("'labels' must be a one or more character strings")
     
  #Compute the widths
    plot.width <- 1 - karyoplot$plot.params$leftmargin - karyoplot$plot.params$rightmargin - length(karyoplot$plot.region)*karyoplot$plot.params$ideogramlateralmargin
    plot.bases <- sum(width(karyoplot$plot.region))
    sw <- graphics::strwidth(s = labels, cex=cex, units = "figure")
    
    sw.in.bases <- (plot.bases/plot.width)*sw
  
  #Compute the heights
    plot.height <- karyoplot$plot.params[[paste0("data", data.panel, "height")]]/karyoplot$plot$ymax #Portion of the total height dedicated to the selected data.panel
    plot.runits <- 1  #tot r units of the data.panel (it's 1 by definition)
    sh <- graphics::strheight(s=labels, cex=cex, units = "figure")
  
    sh.in.runits <- (plot.runits/plot.height)*sh
    
  return(list(width=sw.in.bases, height=sh.in.runits))
}





#change NA's to "" in a char vector
naToEmptyChar <- function(x) {
  x[is.na(x)] <- ""
  return(x)
}
