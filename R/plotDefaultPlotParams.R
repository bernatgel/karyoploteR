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
#' @usage plotDefaultPlotParams(plot.type=2, ...)
#' 
#' @param plot.type (numeric) plot the params of this plot type. Currently, only plot type 2 is accepted. (defaults to 2)
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


plotDefaultPlotParams <- function(plot.type=2, ...) {
  
  valid.plot.types <- c(2) #c(1:4)
  
  if(!plot.type %in% valid.plot.types) {
    stop(paste0("plot.type is not valid. Select a valid value: ", paste0(valid.plot.types, collapse=", ")))
  }
  
  #TODO Implement the plot params for the other plot.types
  if(plot.type == 2) { #Horizontal. Data above and below the ideogram
    
    
    kp <- plotKaryotype(genome="hg19", chromosomes=c("chr21", "chr22"), plot.type=plot.type, ...)
    pp <- kp$plot.params
    
    chrlen <- end(kp$genome[1])
    
    annotateDataPanel <- function(dp, pp, ...) {
        dp.mid = (pp[[paste0("data", dp, "max")]] - pp[[paste0("data", dp, "min")]])/2
        kpDataBackground(kp, data.panel=dp, ...)
        kpText(kp, chr=c("chr21", "chr22"), x=chrlen/2, y=dp.mid, labels = paste0("data.panel=", dp), data.panel=dp, ...)
        
        kpSegments(kp, chr=c("chr21", "chr22"), x0=chrlen/10, x1=(chrlen/10)*1.2, y0=pp[[paste0("data", dp, "max")]], y1=pp[[paste0("data", dp, "max")]], data.panel=dp, ...)  
        kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*1.25, y=pp[[paste0("data", dp, "max")]], labels=paste0("data", dp, "max=", pp[[paste0("data", dp, "max")]]), data.panel=dp, pos=4, ...)  
        
        kpSegments(kp, chr=c("chr21", "chr22"), x0=chrlen/10, x1=(chrlen/10)*1.2, y0=pp[[paste0("data", dp, "min")]], y1=pp[[paste0("data", dp, "min")]], data.panel=dp, ...)  
        kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*1.25, y=pp[[paste0("data", dp, "min")]], labels=paste0("data", dp, "min=", pp[[paste0("data", dp, "min")]]), data.panel=dp, pos=4, ...)  
        
        kpArrows(kp, chr=c("chr21", "chr22"), x0=chrlen/10*8, x1=chrlen/10*8, y0=pp[[paste0("data", dp, "min")]], y1=pp[[paste0("data", dp, "max")]], data.panel=dp, code=3, length=0.07, ...)  
        kpText(kp, chr=c("chr21", "chr22"), x=(chrlen/10)*7.95, y=dp.mid, labels=paste0("data", dp, "height=", pp[[paste0("data", dp, "height")]]), data.panel=dp, pos=2, ...)  
    }
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
      marg.y <- kp$plot$ymax /10 * 9.4
      graphics::arrows(x0 = 0, x1=kp$plot.params$leftmargin, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
      graphics::text(x=kp$plot.params$leftmargin/2, y=marg.y, labels=paste0("leftmargin=", kp$plot.params$leftmargin), pos=3, ...)
      
      graphics::arrows(x0 = kp$plot$xmax - kp$plot.params$rightmargin, x1=kp$plot$xmax, y0=marg.y, y1=marg.y, code=3, length=0.07, ...)
      graphics::text(x= kp$plot$xmax - kp$plot.params$rightmargin/2, y=marg.y, labels=paste0("rightmargin=", kp$plot.params$rightmargin), pos=3, ...)
      
      #TODO: in margins
      
      #TODO: ideogram
    kp$endKpPlot()
        
  } 
  invisible(kp)
}