

#' @export kpBars


kpBars <- function(karyoplot, data=NULL, chr=NULL, x0=NULL, x1=x0, y1=NULL, y0=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  if(!is(karyoplot, "KaryoPlot")) stop(paste0("In kpBars: 'karyoplot' must be a valid 'KaryoPlot' object"))
  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())
  
  
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
  
  
  pp <- prepareParameters4("kpRect", karyoplot=karyoplot, data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  

  x0plot <- ccf(chr=pp$chr, x=pp$x0, data.panel=data.panel)$x
  x1plot <- ccf(chr=pp$chr, x=pp$x1, data.panel=data.panel)$x
  y0plot <- ccf(chr=pp$chr, y=pp$y0, data.panel=data.panel)$y
  y1plot <- ccf(chr=pp$chr, y=pp$y1, data.panel=data.panel)$y
  
  
  rect(xleft=x0plot, xright=x1plot, ytop=y1plot, ybottom=y0plot, ...)
  
}
