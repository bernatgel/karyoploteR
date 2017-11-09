#' kpPlotTranscripts
#' 
#' @description 
#' 
#' Plot gene transcripts on the genome, with options to add strand markers and
#' to differentiate between coding and non-coding exons.
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of ç
#'  karyoploteR. It takes a list with the transcripts, coding and non-coding 
#'  exons and creates a traditional boxes and line representation of the
#'  transcripts. Optionally it can add little arrows along the introns to 
#'  show the transcript strand.
#' 
#' @usage kpPlotTranscripts(karyoplot, data, data.panel=1, r0=NULL, r1=NULL, col="black", border=NULL, avoid.overlapping=TRUE, num.layers=NULL, layer.margin=0.05, clipping=TRUE, ...)
#' 
#' @param karyoplot    (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data    (a \code{GRanges}) A GRanges object with the regions to plot.
#' @param data.panel    (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param r0    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param r1    (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)
#' @param col    (color) The background color of the regions. (defaults to black)
#' @param border    (color) The color used to draw the border of the regions. If NULL, no border is drawn. (defaults to NULL)
#' @param avoid.overlapping    (boolean) Whether overlapping regions should be drawn as stacks (TRUE) on drawing one occluding the other in a single layer (FALSE). (defaults to TRUE)
#' @param num.layers    (numeric) The number of layers the plotting space should be divided into to allow for plotting overlapping regions. The lotting region will be divided into this many pieces regardless if any overlapping regions actually exist. If NULL, the maximum number of regions overlapping a single point in the genome. (defaults to NULL)
#' @param layer.margin    (numeric) The blank space left between layers of regions. (defaults to 0.05)
#' @param clipping  (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...    The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
#'  
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'  
#' @seealso \code{\link{kpPlotGenes}}
#' 
#' @examples
#'  
#'  
#'  
#'  
#'@export kpPlotTranscripts


transcripts <- c(toGRanges("chr1", 100, 1000),
                 toGRanges("chr1", 1500, 3000))
names(transcripts) <- c("T1", "T2")
strand(transcripts) <- c("+", "-")

coding.exons <- list("T1"=c(toGRanges("chr1", 200, 300),
                          toGRanges("chr1", 500, 800),
                          toGRanges("chr1", 900, 950)),
                     "T2"=c(toGRanges("chr1", 2200, 2300),
                          toGRanges("chr1", 2500, 2510),
                          toGRanges("chr1", 2700, 2800)))


non.coding.exons <- list("T1"=c(toGRanges("chr1", 100, 199),
                              toGRanges("chr1", 951, 1000)),
                     "T2"=c(toGRanges("chr1", 1500, 1700),
                          toGRanges("chr1", 1900, 1950),
                          toGRanges("chr1", 2100, 2199),
                          toGRanges("chr1", 2801, 3000)))


data <- list(transcripts=transcripts, coding.exons=coding.exons, non.coding.exons=non.coding.exons)

karyoplot <- plotKaryotype(zoom=toGRanges("chr1", 0, 3200))
kpAddBaseNumbers(karyoplot, tick.dist = 400)

kpPlotTranscripts <- function(karyoplot, data, y0=NULL, y1=NULL, non.coding.exons.height=0.5, 
                              detail.level=2,
                              add.strand.marks=TRUE, mark.height=NULL, mark.width=NULL, mark.distance=5,
                              add.transcript.names=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=1,
                              col="black", coding.exons.col=NULL, coding.exons.border.col=NULL, 
                              non.coding.exons.col=NULL, non.coding.exons.border.col=NULL, 
                              introns.col=NULL, marks.col=NULL, transcript.name.col=NULL,
                              ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, data.panel=1, clipping=TRUE, ...)
  

  #karyoplot
    if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(missing(data)) stop("The parameter 'data' is required")
    #TODO: Check data is valid
  #TODO: Check detail.level is 1 or 2
     
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  #if null, get y0 and y1 to occupy the whole data panel
  if(is.null(y0)) y0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(y1)) y1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]


  #Set the colors if needed
  if(is.null(coding.exons.col)) coding.exons.col <- col
  if(is.null(coding.exons.border.col)) coding.exons.border.col <- col
  if(is.null(non.coding.exons.col)) non.coding.exons.col <- col
  if(is.null(non.coding.exons.col)) non.coding.exons.border.col <- col
  if(is.null(introns.col)) introns.col <- col
  if(is.null(marks.col)) marks.col <- col
  if(is.null(transcript.name.col)) transcript.name.col <- col
  
  #Extend y0 and y1 to the length of the transcripts
  y0 <- rep(y0, length.out=length(data$transcripts))
  y1 <- rep(y1, length.out=length(data$transcripts))


  #If needed, compute the sizes of the strand marks
  if(add.strand.marks) {
    #if params are null, make the strand marks with 1/4 of the transcript height
    if(is.null(mark.height)) {
      mark.height <- 0.25
    }
    if(is.null(mark.width)) {
      #Try to make the marks roughly the same witdh as height IN THE FINAL PLOT irrespective of aspect ration!
      #we need to compute the height of the marker in plot coordinates and
      #then find out the number of bases needed to get the same distance 
      #in the x axis
      #TODO: Could it work with the computed aspect ratio as in https://stat.ethz.ch/pipermail/r-help/2005-October/080598.html
      
      pp <- prepareParameters2(karyoplot = karyoplot, chr=karyoplot$chromosomes[1], data = NULL, x =0, y = c(0, mark.height), ymin = 0, ymax=1, r0=r0, r1=r1, data.panel=data.panel)
      mark.height.in.plot <- karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], y=pp$y[2])$y -
        karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], y=pp$y[1])$y
      plot.region.in.plot <- karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], x=end(karyoplot$plot.region))$x - 
        karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], x=start(karyoplot$plot.region))$x
      total.height <- karyoplot$plot$ymax - karyoplot$plot$ymin
      total.width <- karyoplot$plot$xmax - karyoplot$plot$xmin  
      
      mark.h <- mark.height.in.plot/total.height
      mark.width <- width(karyoplot$plot.region)*(mark.h * total.width)/2
    }
  }

  
  #And start plotting
  for(nt in seq_len(length(data$transcripts))) {
    transcript <- data$transcripts[nt]
    t.y0 <- y0[nt]
    t.y1 <- y1[nt]
    mid.transcript <- t.y0 + (t.y1-t.y0)/2
    strand <- as.character(strand(transcript)[1])
  
    if(detail.level==1) { #plot only boxes
      kpRect(karyoplot, data=transcript, y0=t.y0, y1=t.y1, col=col, border=col, ymax=ymax, ymin=ymin, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping)
      
      if(add.strand.marks) {
        if(width(transcript)>2*mark.width) {
          
          num.marks <- max((width(transcript) - mark.width)%/%(mark.width*mark.distance), 0)
          mark.starts <- start(transcript)+mark.width+mark.width*mark.distance*(c(0, seq_len(num.marks-1)))   
          if(strand!="-") {
            kpSegments(karyoplot, chr=seqnames(transcript), x0=mark.starts+mark.width, x1=mark.starts, y0=mid.transcript, y1=mid.transcript+(mark.height/2), col=marks.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
            kpSegments(karyoplot, chr=seqnames(transcript), x0=mark.starts+mark.width, x1=mark.starts, y0=mid.transcript, y1=mid.transcript-(mark.height/2), col=marks.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
          } else {
            kpSegments(karyoplot, chr=seqnames(transcript), x0=mark.starts, x1=mark.starts+mark.width, y0=mid.transcript, y1=mid.transcript+(mark.height/2), col=marks.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
            kpSegments(karyoplot, chr=seqnames(transcript), x0=mark.starts, x1=mark.starts+mark.width, y0=mid.transcript, y1=mid.transcript-(mark.height/2), col=marks.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
          }
          #NOTE: If transcripts have strand "*", should we just NOT plot any mark even if asked?
        }
      }
    } elsif(detail.level==2) {
      coding.exons <- data$coding.exons[[names(transcript)]]
      non.coding.exons <- data$non.coding.exons[[names(transcript)]]
      
      #coding exons
      kpRect(karyoplot, data=coding.exons, y0=t.y0, y1=t.y1, col=coding.exons.col, border=coding.exons.border.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping)
      
      #non-coding exons 
      nc.y0 <- mid.transcript-(t.y1-t.y0)*non.coding.exons.height/2
      nc.y1 <- mid.transcript+(t.y1-t.y0)*non.coding.exons.height/2
      
      kpRect(karyoplot, data=non.coding.exons, y0=nc.y0, y1=nc.y1, col=non.coding.exons.col, border=non.coding.exons.border.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping)
      
      #introns
      introns <- subtractRegions(transcript, c(coding.exons, non.coding.exons, ignore.mcols=TRUE))
      kpSegments(karyoplot, data=introns, y0=mid.transcript, y1=mid.transcript, col=introns.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, clipping=clipping)
                 
      
      
      
        
    }
    
  }
      
  
   
    if(strand.marks) {
      #if params are null, make the strand marks with 1/4 of the transcript height, width similar or equal to the height, separation, 4 heights. Never overlapping an exons even partially.
      if(is.null(mark.height)) {
        mark.height <- 0.25
      }
      if(is.null(mark.width)) {
        #we need to compute the height of the marker in plot coordinates and
        #then find out the number of bases needed to get the same distance 
        #in the x axis
        pp <- prepareParameters2(karyoplot = karyoplot, chr=karyoplot$chromosomes[1], data = NULL, x =0, y = c(0, mark.height), ymin = 0, ymax=1, r0=r0, r1=r1, data.panel=data.panel)
        mark.height.in.plot <- karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], y=pp$y[2])$y -
          karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], y=pp$y[1])$y
        plot.region.in.plot <- karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], x=end(karyoplot$plot.region))$x - 
          karyoplot$coord.change.function(chr=karyoplot$chromosomes[1], x=start(karyoplot$plot.region))$x
        total.height <- karyoplot$plot$ymax - karyoplot$plot$ymin
        total.width <- karyoplot$plot$xmax - karyoplot$plot$xmin  
        
        mark.h <- mark.height.in.plot/total.height
        mark.width <- width(karyoplot$plot.region)*(mark.h * total.width)/2
      }
      
      #now, find the `positions for the mark starts
      introns.with.mark <- introns[width(introns)>2*mark.width]
      marks.per.intron <- (width(introns.with.mark)-mark.width)%/%(mark.width*mark.distance)
      for(i in seq_len(length(introns.with.mark))) {
        mark.starts <- start(introns.with.mark[i])+mark.width+mark.width*mark.distance*(c(0, seq_len(marks.per.intron[i])))   
        strand <- as.character(strand(transcript)[1])
        if(strand=="+") {
          kpSegments(karyoplot, chr=seqnames(introns.with.mark[i]), x0=mark.starts+mark.width, x1=mark.starts, y0=0.5, y1=0.5+(mark.height/2), col=marks.col, r0=r0, r1=r1, data.panel=data.panel)
          kpSegments(karyoplot, chr=seqnames(introns.with.mark[i]), x0=mark.starts+mark.width, x1=mark.starts, y0=0.5, y1=0.5-(mark.height/2), col=marks.col, r0=r0, r1=r1, data.panel=data.panel)
        } else {
          kpSegments(karyoplot, chr=seqnames(introns.with.mark[i]), x0=mark.starts, x1=mark.starts+mark.width, y0=0.5, y1=0.5+(mark.height/2), col=marks.col, r0=r0, r1=r1, data.panel=data.panel)
          kpSegments(karyoplot, chr=seqnames(introns.with.mark[i]), x0=mark.starts, x1=mark.starts+mark.width, y0=0.5, y1=0.5-(mark.height/2), col=marks.col, r0=r0, r1=r1, data.panel=data.panel)
        }
      }
    }
    
    #plot the transcript name
    if(!is.null(transcript.name)) {
      kpPlotNames()
    }
  }    
  
  
  
  
  
  
  
  invisible(karyoplot)
}


