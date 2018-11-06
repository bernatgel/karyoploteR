#' plotDefaultPlotParameters
#' 
#' @description 
#' 
#' Creates a karyoplot with the default parameters drawn.
#' 
#' @details 
#'  
#'  Given a plot.type, this function creates a new karyoplot with lines and arrows showing 
#'  the meaning and values of the plot.params
#'  
#' @usage plotDefaultPlotParams(plot.type=2, plot.params=NULL, ...)
#' 
#' @param plot.type (numeric) plot the params of this plot type. Currently, only plot types 2 and 3 are accepted. (defaults to 2)
#' @param plot.params (a plot params object) a plot params object such the one returned by \code{\link{getDefaultPlotParams}}. If specified, it will be used to create the plots.
#' @param ... The ellipsis operator can be used to pass any additional graphics parameter
#' 
#' @return
#' 
#' Returns the original karyoplot object, unchanged.
#'    
#' @seealso \code{\link{plotKaryotype}}
#' 
#' 
#' @examples
#' 
#' kp <- plotDefaultPlotParams(plot.type=2)
#'
#'  
#' @export plotDefaultPlotParams
#' 


plotDefaultPlotParams <- function(plot.type=2, plot.params=NULL,  ...) {
  
  #TODO: Implement plot for 4
  
  valid.plot.types <- c(1,2,3,5) #c(1:4)
  
  if(!plot.type %in% valid.plot.types) {
    stop(paste0("plot.type is not valid. Select a valid value: ", paste0(valid.plot.types, collapse=", ")))
  }
  
  annotateDataPanel <- function(dp, pp, ...) {
    dp.mid = (pp[[paste0("data", dp, "max")]] - pp[[paste0("data", dp, "min")]])/2
    dp.03 = (pp[[paste0("data", dp, "max")]] - pp[[paste0("data", dp, "min")]])*0.3
    kpDataBackground(kp, data.panel=dp, ...)
    kpText(kp, chr=c("chr21", "chr22"), x=chrlen/2, y=dp.mid, labels = paste0("data.panel=", dp), data.panel=dp, ...)
    
    kpSegments(kp, chr=c("chr21", "chr22"), x0=chrlen/10, x1=(chrlen/10)*1.2, y0=pp[[paste0("data", dp, "max")]], y1=pp[[paste0("data", dp, "max")]], data.panel=dp, ...)  
    kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*1.25, y=pp[[paste0("data", dp, "max")]], labels=paste0("data", dp, "max=", pp[[paste0("data", dp, "max")]]), data.panel=dp, pos=4, ...)  
    
    kpSegments(kp, chr=c("chr21", "chr22"), x0=chrlen/10, x1=(chrlen/10)*1.2, y0=pp[[paste0("data", dp, "min")]], y1=pp[[paste0("data", dp, "min")]], data.panel=dp, ...)  
    kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*1.25, y=pp[[paste0("data", dp, "min")]], labels=paste0("data", dp, "min=", pp[[paste0("data", dp, "min")]]), data.panel=dp, pos=4, ...)  
    
    kpArrows(kp, chr=c("chr21", "chr22"), x0=chrlen/10*8, x1=chrlen/10*8, y0=pp[[paste0("data", dp, "min")]], y1=pp[[paste0("data", dp, "max")]], data.panel=dp, code=3, length=0.07, ...)  
    kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*7.95, y=dp.03, labels=paste0("data", dp, "height=", pp[[paste0("data", dp, "height")]]), data.panel=dp, pos=2, ...)  
  }
  
  
  #TODO Implement the plot params for the other plot.types
  if(plot.type == 1) { #Horizontal. Data above the ideogram
    if(is.null(plot.params)) {
      plot.params <- getDefaultPlotParams(plot.type=plot.type)
    } 
    
    kp <- plotKaryotype(genome="hg19", chromosomes=c("chr21", "chr22"), plot.type=plot.type, plot.params = plot.params, ...)
    pp <- kp$plot.params
    
    chrlen <- end(kp$genome[1])
    
    
    annotateDataPanel(1, pp, ...)
    
    #mainmargins
    kp$beginKpPlot()
    graphics::abline(h=c(kp$plot.params$bottommargin, kp$plot.params$bottommargin + kp$chromosome.height, kp$plot.params$bottommargin + 2*kp$chromosome.height), ...)
    marg.x <- kp$plot$xmax /10 * 8
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=0, y1=kp$plot.params$bottommargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot.params$bottommargin/2,labels=paste0("bottommargin=", kp$plot.param$bottommargin), pos=2, ...)
    
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=kp$plot$ymax, y1=kp$plot$ymax - kp$plot.params$topmargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot$ymax - kp$plot.params$topmargin/2, labels=paste0("topmargin=", kp$plot.param$topmargin), pos=2, ...)
    
    #chromosome height
    ch.x <- kp$plot$xmax /10 * 9.2
    graphics::arrows(x0 = ch.x, x1=ch.x, y0=kp$plot.params$bottommargin + kp$chromosome.height, y1=kp$plot.params$bottommargin + 2*kp$chromosome.height, code=3, length=0.07, ...)
    graphics::text(x=ch.x, y=kp$plot.params$bottommargin + 1.5*kp$chromosome.height, labels=paste0("chromosome.height (computed)=", kp$chromosome.height), pos=2, ...)
    
    graphics::abline(v=c(kp$plot.params$leftmargin, kp$plot$xmax-kp$plot.params$rightmargin), lty=2, ...)
    marg.y <- kp$plot$ymax /10 * 8.5
    graphics::arrows(x0 = 0, x1=kp$plot.params$leftmargin, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x=kp$plot.params$leftmargin/2, y=marg.y, labels=paste0("leftmargin=", kp$plot.params$leftmargin), pos=3, ...)
    
    graphics::arrows(x0 = kp$plot$xmax - kp$plot.params$rightmargin, x1=kp$plot$xmax, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x= kp$plot$xmax - kp$plot.params$rightmargin/2, y=marg.y, labels=paste0("rightmargin=", kp$plot.params$rightmargin), pos=3, ...)
    
    ideoh.x <- kp$plot$xmax /10 * 5.5
    graphics::arrows(x0 = ideoh.x, x1=ideoh.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideoh.x, y=kp$ideogram.mid("chr21"), labels=paste0("ideogramheight=", pp$ideogramheight), pos=4, ...)
    
    #in and out margins
    ideom.x <- kp$plot$xmax /10 * 7.2
    graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin/2, labels=paste0("data1inmargin=", pp$data1inmargin), pos=4, ...)

    outm.x <- kp$plot$xmax /10 * 7
    outm.y0 <- kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin + pp$data1height
    outm.y1 <- outm.y0 + pp$data1outmargin
    if(pp$data1outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
    graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data1outmargin=", pp$data1outmargin), pos=2, ...)

    kp$endKpPlot()
    
  } 
  if(plot.type == 2) { #Horizontal. Data above and below the ideogram
    
    if(is.null(plot.params)) {
      plot.params <- getDefaultPlotParams(plot.type=2)
    } 
    
    kp <- plotKaryotype(genome="hg19", chromosomes=c("chr21", "chr22"), plot.type=plot.type, plot.params = plot.params, ...)
    pp <- kp$plot.params
    
    chrlen <- end(kp$genome[1])
    
   
    annotateDataPanel(1, pp, ...)
    annotateDataPanel(2, pp, ...)
    
    
    #outmargins
    kp$beginKpPlot()
      graphics::abline(h=c(kp$plot.params$bottommargin, kp$plot.params$bottommargin + kp$chromosome.height, kp$plot.params$bottommargin + 2*kp$chromosome.height), ...)
      marg.x <- kp$plot$xmax /10 * 8
      graphics::arrows(x0 = marg.x, x1=marg.x, y0=0, y1=kp$plot.params$bottommargin, code=3, length=0.07, ...)
      graphics::text(x=marg.x, y=kp$plot.params$bottommargin/2,labels=paste0("bottommargin=", kp$plot.param$bottommargin), pos=2, ...)
      
      graphics::arrows(x0 = marg.x, x1=marg.x, y0=kp$plot$ymax, y1=kp$plot$ymax - kp$plot.params$topmargin, code=3, length=0.07, ...)
      graphics::text(x=marg.x, y=kp$plot$ymax - kp$plot.params$topmargin/2, labels=paste0("topmargin=", kp$plot.param$topmargin), pos=2, ...)
      
      ch.x <- kp$plot$xmax /10 * 9.2
      graphics::arrows(x0 = ch.x, x1=ch.x, y0=kp$plot.params$bottommargin + kp$chromosome.height, y1=kp$plot.params$bottommargin + 2*kp$chromosome.height, code=3, length=0.07, ...)
      graphics::text(x=ch.x, y=kp$plot.params$bottommargin + 1.5*kp$chromosome.height, labels=paste0("chromosome.height (computed)=", kp$chromosome.height), pos=2, ...)
      
      graphics::abline(v=c(kp$plot.params$leftmargin, kp$plot$xmax-kp$plot.params$rightmargin), lty=2, ...)
      marg.y <- kp$plot$ymax /10 * 8.5
      graphics::arrows(x0 = 0, x1=kp$plot.params$leftmargin, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
      graphics::text(x=kp$plot.params$leftmargin/2, y=marg.y, labels=paste0("leftmargin=", kp$plot.params$leftmargin), pos=3, ...)
      
      graphics::arrows(x0 = kp$plot$xmax - kp$plot.params$rightmargin, x1=kp$plot$xmax, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
      graphics::text(x= kp$plot$xmax - kp$plot.params$rightmargin/2, y=marg.y, labels=paste0("rightmargin=", kp$plot.params$rightmargin), pos=3, ...)
      
      ideoh.x <- kp$plot$xmax /10 * 5.5
      graphics::arrows(x0 = ideoh.x, x1=ideoh.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
      graphics::text(x= ideoh.x, y=kp$ideogram.mid("chr21"), labels=paste0("ideogramheight=", pp$ideogramheight), pos=4, ...)
      
      #in margins
      ideom.x <- kp$plot$xmax /10 * 7.2
      graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin, y1=kp$ideogram.mid("chr21") - pp$ideogramheight/2, code=3, length=0.07, ...)
      graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin/2, labels=paste0("data2inmargin=", pp$data2inmargin), pos=4, ...)
      
      graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
      graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin/2, labels=paste0("data1inmargin=", pp$data1inmargin), pos=4, ...)
      
      outm.x <- kp$plot$xmax /10 * 7
      outm.y0 <- kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin + pp$data1height
      outm.y1 <- outm.y0 + pp$data1outmargin
      if(pp$data1outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
      graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data1outmargin=", pp$data1outmargin), pos=2, ...)
      
      outm.x <- kp$plot$xmax /10 * 7
      outm.y0 <- kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin - pp$data2height
      outm.y1 <- ifelse(pp$data2outmargin!=0, outm.y0 - pp$data2outmargin, outm.y0 - pp$data2outmargin - 0.001)
      if(pp$data2outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
      graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data2outmargin=", pp$data2outmargin), pos=2, ...)
      
      
      
      
      kp$endKpPlot()
        
  } 
  
  if(plot.type == 3) { #Horizontal in a single line. Data above and below the ideogram
    
    if(is.null(plot.params)) {
      plot.params <- getDefaultPlotParams(plot.type=3)
      #plot.params$ideogramlateralmargin=0.01
    } 
    
    kp <- plotKaryotype(genome="hg19", chromosomes=c("chr21", "chr22"), plot.type=plot.type, plot.params = plot.params, ...)
    pp <- kp$plot.params

    chrlen <- end(kp$genome[1])
    
    annotateDataPanel(1, pp, ...)
    annotateDataPanel(2, pp, ...)
    
    #outmargins
    kp$beginKpPlot()
    graphics::abline(h=c(kp$plot.params$bottommargin, kp$plot.params$bottommargin + kp$chromosome.height, kp$plot.params$bottommargin + 2*kp$chromosome.height), ...)
    marg.x <- kp$plot$xmax /10 * 8
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=0, y1=kp$plot.params$bottommargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot.params$bottommargin/2,labels=paste0("bottommargin=", kp$plot.param$bottommargin), pos=2, ...)
    
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=kp$plot$ymax, y1=kp$plot$ymax - kp$plot.params$topmargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot$ymax - kp$plot.params$topmargin/2, labels=paste0("topmargin=", kp$plot.param$topmargin), pos=2, ...)
    
    ch.x <- kp$plot$xmax /10 * 9.2
    graphics::arrows(x0 = ch.x, x1=ch.x, y0=kp$plot.params$bottommargin, y1=kp$plot.params$bottommargin + 1*kp$chromosome.height, code=3, length=0.07, ...)
    graphics::text(x=ch.x, y=kp$plot.params$bottommargin + 0.7*kp$chromosome.height, labels=paste0("chromosome.height (computed)=", kp$chromosome.height), pos=2, ...)
    
    graphics::abline(v=c(kp$plot.params$leftmargin, kp$plot$xmax-kp$plot.params$rightmargin), lty=2, ...)
    marg.y <- kp$plot$ymax /10 * 8.5
    graphics::arrows(x0 = 0, x1=kp$plot.params$leftmargin, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x=kp$plot.params$leftmargin/2, y=marg.y, labels=paste0("leftmargin=", kp$plot.params$leftmargin), pos=3, ...)
    
    graphics::arrows(x0 = kp$plot$xmax - kp$plot.params$rightmargin, x1=kp$plot$xmax, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x= kp$plot$xmax - kp$plot.params$rightmargin/2, y=marg.y, labels=paste0("rightmargin=", kp$plot.params$rightmargin), pos=3, ...)
    
    #ideogram
    ideoh.x <- kp$plot$xmax /10 * 5.5
    graphics::arrows(x0 = ideoh.x, x1=ideoh.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideoh.x, y=kp$ideogram.mid("chr21"), labels=paste0("ideogramheight=", pp$ideogramheight), pos=4, ...)
    
    #in margins
    ideom.x <- kp$plot$xmax /10 * 7.2
    graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin, y1=kp$ideogram.mid("chr21") - pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin/2, labels=paste0("data2inmargin=", pp$data2inmargin), pos=4, ...)
    
    graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin/2, labels=paste0("data1inmargin=", pp$data1inmargin), pos=4, ...)
    
    outm.x <- kp$plot$xmax /10 * 7
    outm.y0 <- kp$ideogram.mid("chr21") + pp$ideogramheight/2 + pp$data1inmargin + pp$data1height
    outm.y1 <- ifelse(pp$data1outmargin!=0, outm.y0 + pp$data1outmargin, outm.y0 + pp$data1outmargin + 0.001)
    if(pp$data1outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
    graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data1outmargin=", pp$data1outmargin), pos=2, ...)
    
    outm.x <- kp$plot$xmax /10 * 7
    outm.y0 <- kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin - pp$data2height
    outm.y1 <- ifelse(pp$data2outmargin!=0, outm.y0 - pp$data2outmargin, outm.y0 - pp$data2outmargin - 0.001)
    if(pp$data2outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
    graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data2outmargin=", pp$data2outmargin), pos=2, ...)
    
    
    #horizontal margin
    hmarg.y <- kp$coord.change.function(chr="chr21", y=pp$data1max*0.8, data.panel = 1)$y
    hmarg.x0 <- kp$coord.change.function(chr = "chr21", x = kp$chromosome.lengths[1], data.panel = 1)$x
    hmarg.x1 <- ifelse(pp$ideogramlateralmargin!=0, hmarg.x0 + pp$ideogramlateralmargin, hmarg.x0 + pp$ideogramlateralmargin+0.001)
    graphics::arrows(x0 = hmarg.x0, x1=hmarg.x1, y0=hmarg.y, y1=hmarg.y, code=3, length=0.03, ...)
    graphics::text(x= mean(hmarg.x0, hmarg.x1), y=hmarg.y, labels=paste0("ideogramlateralmargin=", pp$ideogramlateralmargin), pos=3, ...)

    
    kp$endKpPlot()
    
  } 
  
  
  
  if(plot.type == 5) { #Horizontal in a single line. Data above and below the ideogram
    
    if(is.null(plot.params)) {
      plot.params <- getDefaultPlotParams(plot.type=5)
      #plot.params$ideogramlateralmargin=0.01
    } 
    
    kp <- plotKaryotype(genome="hg19", chromosomes=c("chr21", "chr22"), plot.type=plot.type, plot.params = plot.params, ...)
    pp <- kp$plot.params
    
    chrlen <- end(kp$genome[1])
    
    annotateDataPanel(2, pp, ...)
    
    #outmargins
    kp$beginKpPlot()
    graphics::abline(h=c(kp$plot.params$bottommargin, kp$plot.params$bottommargin + kp$chromosome.height), ...)
    marg.x <- kp$plot$xmax /10 * 8
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=0, y1=kp$plot.params$bottommargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot.params$bottommargin/2,labels=paste0("bottommargin=", kp$plot.param$bottommargin), pos=2, ...)
    
    graphics::arrows(x0 = marg.x, x1=marg.x, y0=kp$plot$ymax, y1=kp$plot$ymax - kp$plot.params$topmargin, code=3, length=0.07, ...)
    graphics::text(x=marg.x, y=kp$plot$ymax - kp$plot.params$topmargin/2, labels=paste0("topmargin=", kp$plot.param$topmargin), pos=2, ...)
    
    ch.x <- kp$plot$xmax /10 * 9.2
    graphics::arrows(x0 = ch.x, x1=ch.x, y0=kp$plot.params$bottommargin, y1=kp$plot.params$bottommargin + 1*kp$chromosome.height, code=3, length=0.07, ...)
    graphics::text(x=ch.x, y=kp$plot.params$bottommargin + 0.7*kp$chromosome.height, labels=paste0("chromosome.height (computed)=", kp$chromosome.height), pos=2, ...)
    
    graphics::abline(v=c(kp$plot.params$leftmargin, kp$plot$xmax-kp$plot.params$rightmargin), lty=2, ...)
    marg.y <- kp$plot$ymax /10 * 8.5
    graphics::arrows(x0 = 0, x1=kp$plot.params$leftmargin, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x=kp$plot.params$leftmargin/2, y=marg.y, labels=paste0("leftmargin=", kp$plot.params$leftmargin), pos=3, ...)
    
    graphics::arrows(x0 = kp$plot$xmax - kp$plot.params$rightmargin, x1=kp$plot$xmax, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
    graphics::text(x= kp$plot$xmax - kp$plot.params$rightmargin/2, y=marg.y, labels=paste0("rightmargin=", kp$plot.params$rightmargin), pos=3, ...)
    
    #ideogram
    ideoh.x <- kp$plot$xmax /10 * 5.5
    graphics::arrows(x0 = ideoh.x, x1=ideoh.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2, y1=kp$ideogram.mid("chr21") + pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideoh.x, y=kp$ideogram.mid("chr21"), labels=paste0("ideogramheight=", pp$ideogramheight), pos=4, ...)
    
    #in margins
    ideom.x <- kp$plot$xmax /10 * 7.2
    graphics::arrows(x0 = ideom.x, x1=ideom.x, y0=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin, y1=kp$ideogram.mid("chr21") - pp$ideogramheight/2, code=3, length=0.07, ...)
    graphics::text(x= ideom.x, y=kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin/2, labels=paste0("data2inmargin=", pp$data2inmargin), pos=4, ...)
    
    outm.x <- kp$plot$xmax /10 * 7
    outm.y0 <- kp$ideogram.mid("chr21") - pp$ideogramheight/2 - pp$data2inmargin - pp$data2height
    outm.y1 <- ifelse(pp$data2outmargin!=0, outm.y0 - pp$data2outmargin, outm.y0 - pp$data2outmargin - 0.001)
    if(pp$data2outmargin!=0) graphics::arrows(x0 = outm.x, x1=outm.x, y0=outm.y0, y1=outm.y1, code=3, length=0.07, ...)
    graphics::text(x= outm.x, y=mean(outm.y0, outm.y1), labels=paste0("data2outmargin=", pp$data2outmargin), pos=2, ...)
    
    
    #horizontal margin
    hmarg.y <- kp$coord.change.function(chr="chr21", y=pp$data2max*0.8, data.panel = 2)$y
    hmarg.x0 <- kp$coord.change.function(chr = "chr21", x = kp$chromosome.lengths[1], data.panel=1)$x
    hmarg.x1 <- ifelse(pp$ideogramlateralmargin!=0, hmarg.x0 + pp$ideogramlateralmargin, hmarg.x0 + pp$ideogramlateralmargin+0.001)
    graphics::arrows(x0 = hmarg.x0, x1=hmarg.x1, y0=hmarg.y, y1=hmarg.y, code=3, length=0.03, ...)
    graphics::text(x= mean(hmarg.x0, hmarg.x1), y=hmarg.y, labels=paste0("ideogramlateralmargin=", pp$ideogramlateralmargin), pos=3, ...)
    
    
    kp$endKpPlot()
    
  } 
  
  
  invisible(kp)
}