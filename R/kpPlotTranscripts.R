#' kpPlotTranscripts
#' 
#' @description 
#' 
#' Plot gene transcripts on the genome, with options to add strand markers and
#' to differentiate between coding and non-coding exons.
#' 
#' @details 
#'  
#'  This is one of the high-level, or specialized, plotting functions of
#'  karyoploteR. It takes a list with the transcripts, coding and non-coding 
#'  exons and creates a traditional boxes and line representation of the
#'  transcripts. With y0 and y1, it is possible to specify a different vertical 
#'  position and different height for each transcript. 
#'  It can add little arrows,
#'  strand marks, along the introns to show the transcript strand. The marks 
#'  appearance can be customized using 4 different parameters specifying the 
#'  height (relative to the height of the transcript), width of each mark 
#'  (relative to the height), distance between marks (relative to the width) and
#'  color of the marks. Marks are centered on the space they have available and
#'  if the available space to too tight for a single mark, no mark will be 
#'  plotted. The direction of the marks is based on the transcript strand as
#'  specified in \code{data$transcripts} object.
#'  Two detail levels are available: \code{detail.level=1} will represent 
#'  transcripts as solid boxes (optionally with strand marks along the whole
#'  transcript); \code{detail.level=2} will represent the internal structure of 
#'  the transcripts -coding and non coding exons, introns, and optionally the
#'  strand marks only in the introns.
#'  
#' @note IMPORTANT: The direction of the strand marks is taken from the strand
#' information in the \code{data$transcripts} object. If transcripts have no
#' strand information there (they have \code{strand="*"}), no marks will be 
#' drawn.
#' 
#' 
#' 
#' @usage kpPlotTranscripts(karyoplot, data, y0=NULL, y1=NULL, non.coding.exons.height=0.5, 
#' detail.level=2,
#' add.strand.marks=TRUE, mark.height=0.20, mark.width=1, mark.distance=4,
#' add.transcript.names=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=1,
#' col="black", border=NULL, coding.exons.col=NULL, coding.exons.border.col=NULL, 
#' non.coding.exons.col=NULL, non.coding.exons.border.col=NULL, 
#' introns.col=NULL, marks.col=NULL, transcript.name.col=NULL,
#' ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, autotrack=NULL, 
#' data.panel=1, clipping=TRUE, ...) 
#' 
#' @param karyoplot (a \code{KaryoPlot} object) This is the first argument to all data plotting functions of \code{karyoploteR}. A KaryoPlot object referring to the currently active plot.
#' @param data  (a \code{list}) data must be a list with an element called 'transcripts' with a \code{GRanges} with all transcripts to be plotted. Additionally, to plot the transcript structure it needs two other elements, 'coding.exons' and 'non.coding.exons', each a list with an element for every transcript with \code{GRanges} objects representing the transcript exons. An example of the data structure needed can be seen in the function examples.
#' @param y0  (numeric) The bottom of the transcripts in the y axis. It can have a different value for each transcript and values will be recycled if needed. If null, it will be set to the minimum y value in the data.panel, usually 0. (Defaults to NULL)
#' @param y1  (numeric) The top of the transcripts in the y axis. It can have a different value for each transcript and values will be recycled if needed. If null, it will be set to the maximum y value in the data.panel, usually 1. (Defaults to NULL)
#' @param non.coding.exons.height  (numeric) The height of the non.coding exons relative to the transcript height. For example, if 0.5, non-coding exons will have a height half the size of the coding ones. (default 0.5) 
#' @param detail.level (numeric: 1 or 2) The detail level of the transcript representation: 1 will plot only boxes representing the transcripts, 2 will plot detailed structure of the trasncripts (coding and non-coding exons and introns). (Defaults to 2)
#' @param add.strand.marks  (boolean) Whether strand marks should be plotted or not. Strand marks are small arrows along the introns (or whole transcripts if detail level=1). (defaults to TRUE)
#' @param mark.height (numeric) The height of the strand marks in "coding exons heights", that is, if mark.height is 0.5, the mark will have a height of half the height of an exon. (defaults to 0.2)
#' @param mark.width  (numeric) The width of the strand marks, in mark heights. mark.width=1 will produce arrow heads with a slope pf 45 degrees. A value higher than 1 will produce smaller angles and a value below 1 larger angles with more vertical lines. (defaults to 1, 45 degrees)
#' @param mark.distance  (numeric) The distance between marks, in mark widths. A distance of 2, will add a space of 2*mark.width between consecutive marks. (defaults to 4)
#' @param add.transcript.names (boolean) Whether to add transcript names to the plot (defailts to TRUE)
#' @param transcript.names  (named character) A named character vector with the labels of the transcripts. If not null, it will be used as a dictionary, so transcript ids should be names and desired labels the values. If NULL, the transcript ids will be used as labels. (defaults to null)
#' @param transcript.name.position  (character) The position of the text relative to the rectangle. Can be "left", "right", "top", "bottom" or "center". (Defaults to "left")
#' @param transcript.name.cex  (numeric) The cex value to plot the transcript labels. (defaults to 1)
#' @param col  (color) The color of the transcripts and labels. It is possible to specify different colors for each element class (transcript names, exons, strand marks...). All elements with no explicit color will be plotted using col. (defaults to "black) 
#' @param border  (color) The color of the border of rectangles representing genes, transcripts and exons. Every element class may have its own specific color using the appropiate parameters. The ones with no explicit color will use border. At the same time, if border is NULL, it will default to col. (fesults to NULL)
#' @param coding.exons.col  (color) The fill color of the rectangles representing the coding exons. If NULL, it will use col. (defaults to NULL)
#' @param coding.exons.border.col  (color) The color of the border of the coding exons. If NULL, it will use border. (defaults to NULL)
#' @param non.coding.exons.col  (color) The fill color of the rectangles representing the non-coding exons. If NULL, it will use col. (defaults to NULL)
#' @param non.coding.exons.border.col  (color) The color of the border of the non-coding exons. If NULL, it will use border. (defaults to NULL)
#' @param introns.col (color) The color of the lines representing the introns. If NULL, it will use col. (defaults to NULL) 
#' @param marks.col   (color) The color of the arrows representing the strand. If NULL, it will use col. (defaults to NULL) 
#' @param transcript.name.col  (color) The color of the transcript labels. If NULL, it will use col. (defaults to NULL)
#' @param ymax  (numeric) The maximum value of \code{y} to be plotted. If NULL, it is set to the max value of the selected data panel. (defaults to NULL)
#' @param ymin  (numeric) The minimum value of \code{y} to be plotted. If NULL, it is set to the min value of the selected data panel. (defaults to NULL) 
#' @param r0  (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL) 
#' @param r1  (numeric) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL) 
#' @param autotrack  (list of numerics) a list numerics with 2 or 3 elements. The first element is the tracks to use with the current plot, the second element is the total number of tracks and the third element is the margin to leave over each track. If the first element, the current track, has more than one element, the plot will span from track `min(autotrack[[1]])` to track `max(autotrack[[1]])`. The margin is specified as the part of a track, by default 0.05, 5 percent of the track height. If NULL, no autotracks will be used. (defaults to NULL)
#' @param data.panel (numeric) The identifier of the data panel where the data is to be plotted. The available data panels depend on the plot type selected in the call to \code{\link{plotKaryotype}}. (defaults to 1)
#' @param clipping (boolean) Only used if zooming is active. If TRUE, the data representation will be not drawn out of the drawing area (i.e. in margins, etc) even if the data overflows the drawing area. If FALSE, the data representation may overflow into the margins of the plot. (defaults to TRUE)
#' @param ...  The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to the R base plotting functions. 
#' 
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
#' #Build the data objet expected by transcripts
#' 
#' transcripts <- c(toGRanges("chr1", 100, 1000),
#'                  toGRanges("chr1", 1500, 3000))
#' names(transcripts) <- c("T1", "T2")
#' strand(transcripts) <- c("+", "-")
#' 
#' coding.exons <- list("T1"=c(toGRanges("chr1", 200, 300),
#'                             toGRanges("chr1", 500, 800),
#'                             toGRanges("chr1", 900, 950)),
#'                      "T2"=c(toGRanges("chr1", 2200, 2300),
#'                             toGRanges("chr1", 2500, 2510),
#'                             toGRanges("chr1", 2700, 2800)))
#' 
#' 
#' non.coding.exons <- list("T1"=c(toGRanges("chr1", 100, 199),
#'                                 toGRanges("chr1", 951, 1000)),
#'                          "T2"=c(toGRanges("chr1", 1500, 1700),
#'                                 toGRanges("chr1", 1900, 1950),
#'                                 toGRanges("chr1", 2100, 2199),
#'                                 toGRanges("chr1", 2801, 3000)))
#' 
#' data <- list(transcripts=transcripts, coding.exons=coding.exons, non.coding.exons=non.coding.exons)
#' 
#' #Create a simple example plot
#' karyoplot <- plotKaryotype(zoom=toGRanges("chr1", 0, 3200))
#' kpAddBaseNumbers(karyoplot, tick.dist = 400)
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0, r1=0.1)
#' 
#' 
#' #Create a plot with different variants of the transcripts
#' karyoplot <- plotKaryotype(zoom=toGRanges("chr1", 0, 3200))
#' kpAddBaseNumbers(karyoplot, tick.dist = 400)
#' #Standard
#' kpPlotTranscripts(karyoplot, data=data, r0=0, r1=0.1)
#' #Customize colors
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.11, r1=0.21, transcript.name.position = "right", non.coding.exons.col = "#888888", non.coding.exons.border.col = "#888888", marks.col="#CC6666", transcript.name.col="#CC6666")
#' #Change vertical position and transcript height
#' kpPlotTranscripts(karyoplot, data=data, y0=c(0, 0.4), y1=c(0.2, 1), r0=0.25, r1=0.8, add.transcript.names = TRUE, transcript.name.position = "left", add.strand.marks = TRUE, clipping=FALSE)
#' #Change detail level, colors and transcript names
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.9, r1=1, detail.level = 1, add.transcript.names = TRUE, transcript.name.position = "top", add.strand.marks = TRUE, clipping=FALSE, col="blue", marks.col="white", transcript.names=list(T1="Transcript 1", T2="Transcript 2"))
#'
#' #Create a plot with different variants of the strand marks
#' karyoplot <- plotKaryotype(zoom=toGRanges("chr1", 0, 3200))
#' kpAddBaseNumbers(karyoplot, tick.dist = 400)
#' #Standard
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0, r1=0.1)
#' #No marks
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.15, r1=0.25, add.strand.marks=FALSE)
#' #Change the strand marks height
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.3, r1=0.4, mark.height=1)
#' #Change the mark width
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.45, r1=0.55, mark.width=2)
#' #Change the mark distance
#' kpPlotTranscripts(karyoplot, data=data, y0=0, y1=1, r0=0.6, r1=0.7, mark.distance=1.5)
#'
#'  
#'  
#' @export kpPlotTranscripts
#'
#' @importFrom graphics par segments




kpPlotTranscripts <- function(karyoplot, data, y0=NULL, y1=NULL, non.coding.exons.height=0.5, 
                              detail.level=2,
                              add.strand.marks=TRUE, mark.height=0.20, mark.width=1, mark.distance=4,
                              add.transcript.names=TRUE, transcript.names=NULL, transcript.name.position="left", transcript.name.cex=1,
                              col="black", border=NULL, coding.exons.col=NULL, coding.exons.border.col=NULL, 
                              non.coding.exons.col=NULL, non.coding.exons.border.col=NULL, 
                              introns.col=NULL, marks.col=NULL, transcript.name.col=NULL,
                              ymax=NULL, ymin=NULL, r0=NULL, r1=NULL, 
                              autotrack=NULL, data.panel=1, clipping=TRUE, ...) {
  

  #karyoplot
    if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
    if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
    if(missing(data)) stop("The parameter 'data' is required")
    #TODO: Check data is valid
  #detail.level
    if(is.null(detail.level)) detail.level <- 2
    if(!(detail.level %in% c(1,2))) warning("Detail level must be 1 (only boxes) or 2 (whole transcript structure). Setting it to 2.")
  
     
  #if null, get the r0 and r1
  if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
  
  #if null, get y0 and y1 to occupy the whole data panel
  if(is.null(y0)) y0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
  if(is.null(y1)) y1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]


  #Set the colors if needed
  if(is.null(border)) border <- col
  if(is.null(coding.exons.col)) coding.exons.col <- col
  if(is.null(coding.exons.border.col)) coding.exons.border.col <- border
  if(is.null(non.coding.exons.col)) non.coding.exons.col <- col
  if(is.null(non.coding.exons.col)) non.coding.exons.border.col <- border
  if(is.null(introns.col)) introns.col <- col
  if(is.null(marks.col)) marks.col <- col
  if(is.null(transcript.name.col)) transcript.name.col <- col
  
  #Extend y0 and y1 to the length of the transcripts
  y0 <- rep(y0, length.out=length(data$transcripts))
  y1 <- rep(y1, length.out=length(data$transcripts))

  #And start plotting
  for(nt in seq_len(length(data$transcripts))) {
    transcript <- data$transcripts[nt]
    chr <- as.character(seqnames(transcript))
    t.y0 <- y0[nt]
    t.y1 <- y1[nt]
    mid.transcript <- t.y0 + (t.y1-t.y0)/2
    transcript.height <- t.y1-t.y0
  
    if(detail.level==1) { #plot only boxes
      kpRect(karyoplot, data=transcript, y0=t.y0, y1=t.y1, col=col, border=border, ymax=ymax, ymin=ymin, r0=r0, r1=r1, autotrack = autotrack, data.panel=data.panel, clipping=clipping)
      
      if(add.strand.marks==TRUE) {
        addStrandMarks(karyoplot, transcript, transcript.height=transcript.height, mid.transcript = mid.transcript, mark.height=mark.height, mark.width=mark.width, mark.distance=mark.distance, marks.col=marks.col, r0=r0, r1=r1, autotrack=autotrack, ymin=ymin, ymax=ymax, data.panel=data.panel, clipping=clipping)
      }
    } else if(detail.level==2) {
      coding.exons <- data$coding.exons[[names(transcript)]]
      non.coding.exons <- data$non.coding.exons[[names(transcript)]]
      
      #Plot introns, non coding exons and coding exons in that order so overlapping margins do not enter the larger figures

      #introns
      introns <- setdiff(transcript, c(coding.exons, non.coding.exons, ignore.mcols=TRUE), ignore.strand=TRUE)
      strand(introns) <- strand(transcript)

      kpSegments(karyoplot, data=introns, y0=mid.transcript, y1=mid.transcript, col=introns.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, autotrack=autotrack, data.panel=data.panel, clipping=clipping)

      if(add.strand.marks==TRUE) {
        addStrandMarks(karyoplot, introns, transcript.height=transcript.height, mid.transcript = mid.transcript, mark.height=mark.height, mark.width=mark.width, mark.distance=mark.distance, marks.col=marks.col, r0=r0, r1=r1, autotrack=autotrack, ymin=ymin, ymax=ymax, data.panel=data.panel, clipping=clipping)
      }
      
            
      #non-coding exons 
      nc.y0 <- mid.transcript-(t.y1-t.y0)*non.coding.exons.height/2
      nc.y1 <- mid.transcript+(t.y1-t.y0)*non.coding.exons.height/2
      
      kpRect(karyoplot, data=non.coding.exons, y0=nc.y0, y1=nc.y1, col=non.coding.exons.col, border=non.coding.exons.border.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, autotrack=autotrack, data.panel=data.panel, clipping=clipping)
      
      #coding exons
      kpRect(karyoplot, data=coding.exons, y0=t.y0, y1=t.y1, col=coding.exons.col, border=coding.exons.border.col, ymin=ymin, ymax=ymax, r0=r0, r1=r1, autotrack=autotrack, data.panel=data.panel, clipping=clipping)
      
    }
    
    #plot the transcript name
    if(add.transcript.names==TRUE) {
      if(!is.null(transcript.names)) {
        transcript.labels <- transcript.names[names(transcript)]
      } else {
        transcript.labels <- names(transcript)
      }
      kpPlotNames(karyoplot, data = transcript, y0 = t.y0, y1=t.y1, labels = transcript.labels, position = transcript.name.position, cex=transcript.name.cex, col=transcript.name.col, r0=r0, r1=r1, autotrack = autotrack, ymin=ymin, ymax=ymax, clipping=clipping, data.panel = data.panel)
    }
  }    
  
  invisible(karyoplot)
}

#Note: the code assumes all regs are in the same chromosome and the same strand
addStrandMarks <- function(karyoplot, regs, transcript.height, mid.transcript,
                           mark.height, mark.width=NULL, mark.distance=NULL, marks.col=NULL,
                           r0=NULL, r1=NULL, autotrack=NULL,
                           ymin=NULL, ymax=NULL, data.panel=NULL, clipping=NULL) {
  #if there are no regions to add marks to, simply return
  if(length(regs)==0) return()
  #since the strand marks are plotted sith direct calls to base R graphics, if clipping is TRUE, cut them with the plotting region
  #intersect removes with ignore.strand==TRUE removes the strand, so get it before
  st <- as.character(strand(regs)[1])
  if(clipping==TRUE) {
    regs <- intersect(regs, karyoplot$plot.region, ignore.strand=TRUE)
    strand(regs) <- st 
  }
  #if after clipping there are no regions to add marks to, simply return
  if(length(regs)==0) return()
  #if the marks would be invisible (height=0, simply return)
  if(mark.height==0) return()
  
  #We'll be plotting something, so begin kpPlot
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  #Assume all regs are in the same chromosome
  chr <- as.character(seqnames(regs))[1]
  
  #Transform mid.transcript and transcript.height to account for r0/r1, ymin, ymax, etc..
  pp <- prepareParameters2(karyoplot, function.name = "addStrandMarks", data = NULL,
                           chr = rep(chr, 3), x=c(0,0,0), y=c(mid.transcript, 0, transcript.height),
                           r0=r0, r1=r1, autotrack=autotrack, ymin=ymin, ymax=ymax, data.panel = data.panel)
  mid.transcript <- pp$y[1]
  transcript.height <- abs(pp$y[2]-pp$y[3])
  
  #This function works in plot coordinates, to allow for aspect ratio corrections
  ccf <- karyoplot$coord.change.function
  
  #strand marks height
  t.mark.height <- transcript.height*mark.height/2
  plot.mark.height <- diff(ccf(chr = c(chr, chr), x=0, y=c(0, t.mark.height))$y)
  
  #Compute the mark.width relative to the mark.height. 
  #Adjust for the aspcet ratio of the plot coordinates (in points)
  asp.usr <- diff(par("usr")[1:2])/diff(par("usr")[3:4])
  #Adjust for the aspect ratio of the canvas (in "inches")
  asp.pin <- par("pin")[2]/par("pin")[1]
  
  plot.mark.width <- plot.mark.height*asp.usr*asp.pin*mark.width
  

  for(nr in seq_len(length(regs))) {
    reg <- regs[nr]
    if(st=="*") {
      #NOTE: If transcripts have strand "*", should we just NOT plot any mark
      return()
    }
    #Now, given the regions, determine if how many marks we should plot and where
    plot.coords <- ccf(chr=c(chr, chr), x=c(start(reg), end(reg)), y=c(mid.transcript, mid.transcript))
    
    reg.width <- diff(plot.coords$x)
    if(reg.width > 1.4*plot.mark.width) {
      num.marks <- (reg.width-plot.mark.width)%/%(plot.mark.width*mark.distance)+1
      left.margin <- (reg.width-((num.marks-1)*plot.mark.width*mark.distance)-plot.mark.width)/2
      mark.starts <- plot.coords$x[1]+left.margin+plot.mark.width/2+plot.mark.width*mark.distance*c(0:(num.marks-1))
      
      if(st=="-") {
        segments(x0 = mark.starts-plot.mark.width/2, x1=mark.starts+plot.mark.width/2, y0=plot.coords$y[1], y1=plot.coords$y[1]+plot.mark.height, col=marks.col)
        segments(x0 = mark.starts-plot.mark.width/2, x1=mark.starts+plot.mark.width/2, y0=plot.coords$y[1], y1=plot.coords$y[1]-plot.mark.height, col=marks.col)
      } else {
        segments(x0 = mark.starts-plot.mark.width/2, x1=mark.starts+plot.mark.width/2, y0=plot.coords$y[1]+plot.mark.height, y1=plot.coords$y[1], col=marks.col)
        segments(x0 = mark.starts-plot.mark.width/2, x1=mark.starts+plot.mark.width/2, y0=plot.coords$y[1]-plot.mark.height, y1=plot.coords$y[1], col=marks.col)
      }
    }
  }
}


