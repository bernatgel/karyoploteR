
# 
# 
# kpPlotWidth <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, col="blue", ...) {
#   
#   if(is.null(r0)) r0 <- karyoplot$plot.params[[paste0("data", data.panel, "min")]]
#   if(is.null(r1)) r1 <- karyoplot$plot.params[[paste0("data", data.panel, "max")]]
#   
#   app <- prepareParameters2("kpLines", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, ymin=ymin, ymax=ymax, r0=(r0+r1)/2, r1=r1, data.panel=data.panel, ...)
#   bpp <- prepareParameters2("kpLines", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, ymin=ymin, ymax=ymax, r0=(r0+r1)/2, r1=r0, data.panel=data.panel, ...)
#   ccf <- karyoplot$coord.change.function
#     
#   kpPlot2Lines(karyoplot=karyoplot, 
#                achr=app$chr, ax=app$x, ay=app$y, acol=col, afill=NULL,
#                bchr=bpp$chr, bx=bpp$x, by=bpp$y, bcol=col, bfill=NULL,                
#                ymin=NULL, ymax=NULL, r0=r0, r1=r1, data.panel=data.panel, ...)
#   
# }
