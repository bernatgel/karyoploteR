#' kpAddColorRect
#' 
#' @description 
#' 
#' Add color rectangles next to the data panels. Ideal to identify the data in the plot
#' 
#' @details 
#' 
#' Given a KaryoPlot object, plot colored rectangles on the side of the data panels to help identify the different samples or types of data plotted
#' 
#' @usage kpAddColorRect(karyoplot, col="gray", rect.width=0.02, rect.margin=0.01,  side="left", y0=0, y1=1, r0=NULL, r1=NULL, data.panel=1, border=NA, ...)
#' 
#' @param karyoplot    a \code{karyoplot} object returned by a call to \code{plotKaryotype}
#' @param col   (color) the color of the rectangle
#' @param rect.width    (numeric) the width of the rectangle in plot coordinates (the whole plot has a width of 1). Usual value might be 0.05. Can be negative. (defaults to 0.02)
#' @param rect.margin    (numeric) the additional the margin between the rectangle and the chromosome. In plot coordinates (the whole plot has a width of 1). Usual value might be 0.05. Can be negative. (defaults to 0.01)
#' @param side ("left" or "right") The side of the plot where to plot the labels. (defaults to "left") 
#' @param y0 (numeric) If the vertical space defined by r0 and r1 goes from 0 to 1, the value at which the bottom of the rectangle starts. Ex: y0=0.5, y1=1 will create a rectangle occupying the top half of the total allocated vertical space. (defaults to 0)
#' @param y1 (numeric) If the vertical space defined by r0 and r1 goes from 0 to 1, the value at which the top of the rectangle ends. Ex: y0=0, y1=0.5 will create a rectangle occupying the bottom half of the total allocated vertical space. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the rectangle. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to position the rectangle. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param data.panel    (numeric) The identifier of the data panel where the labels are to be added. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param border (color) The color of the rectangle border. If NA, no border will be plotted. (defaults to NA)
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#' 
#' @return
#' invisibly returns the given karyoplot object
#'  
#' @seealso \code{\link{plotKaryotype}}, \code{\link{autotrack}}
#' 
#' @examples
#'
#'  num.samples <- 10
#' samples.metadata <- data.frame(stage=c("I", "I", "I", "III", "IV", "III", "II", "IV", "IV", "IV"),
#'                                 metastasis=c(0,0,0,0,1,0,0,0,1,0),
#'                                 subtype=c("A", "A", "B", "B", "A", "B", "B", "B", "A", "B"))
#' 
#' kp <- plotKaryotype(plot.type=4)
#' for(i in 1:num.samples) {
#'   #metastasis
#'   kpAddColorRect(kp, col = colByCategory(samples.metadata$metastasis, color=c("white", "black"))[i],
#'                  rect.margin = 0.01, rect.width = 0.01, r0=autotrack(i, num.samples))
#'   #stage
#'   kpAddColorRect(kp, 
#'     col = colByCategory(samples.metadata$stage, colors=c("I"="plum1", "II"="hotpink", "III"="orange", "IV"="red2"))[i],
#'     rect.margin = 0.021, rect.width = 0.01, r0=autotrack(i, num.samples))
#'   #subtype
#'   kpAddColorRect(kp, 
#'                  col = colByCategory(samples.metadata$subtype, colors=c("dodgerblue", "gold"))[i],
#'                  rect.margin = 0.032, rect.width = 0.01, r0=autotrack(i, num.samples))
#' }
#' 
#' @export kpAddColorRect
#' 

kpAddColorRect <- function(karyoplot, col="gray", rect.width=0.02, rect.margin=0.01,  side="left", y0=0, y1=1, r0=NULL, r1=NULL, data.panel=1, border=NA, ...) {
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #Depending on the plot type, rectangles will be added to every chromosome or just to the first/last chromosome
  if(karyoplot$plot.type %in% c(1,2,6)) {
    chrs <- karyoplot$chromosomes
  } 
  if(karyoplot$plot.type %in% c(3,4,5,7)) {
    if(side=="left") {
      chrs <- karyoplot$chromosomes[1]
    } else {
      chrs <- karyoplot$chromosomes[length(karyoplot$chromosomes)]
    }
  }
  
  #Compute the vertical positioning
    # 1 - adjust the y in chromosomal coordinates to take r's, y values, etc into account
      adj.y <- prepareParameters2("kpAddColorRect", karyoplot, data=NULL, chr=chrs, 
                    x=0, y=c(y0,y1), ymax=1, ymin=0, r0=r0, r1=r1, data.panel=data.panel)$y
    # 2 - Transform the y values into plot coordinates
      y0.plot <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[1], length(chrs)), data.panel = data.panel)$y
      y1.plot <- karyoplot$coord.change.function(chr = chrs, x = 0, y=rep(adj.y[2], length(chrs)), data.panel = data.panel)$y
  
  #Now compute the horizontal positioning
    if(side == "left") { #if side="left", plot them on the left
      #The right end of the left-side margin minus the rect.margin
      x1.plot <- rep(karyoplot$plot.params$leftmargin, length(chrs)) - rect.margin
      x0.plot <- x1.plot - rect.width
    } else {  #if side="right", plot them on the right
      #The last position of each chrs plus the label.margin
      x0.plot <- karyoplot$coord.change.function(chr = chrs, 
                                           x = end(karyoplot$genome[chrs]),
                                           y=rep(adj.y[1], length(chrs)), 
                                           data.panel = data.panel)$x + rect.margin
      x1.plot <- x0.plot + rect.width
    }

  graphics::rect(xleft = x0.plot, xright = x1.plot, ybottom = y0.plot, ytop = y1.plot, col = col, border=border, ...)

  invisible(karyoplot)
}

