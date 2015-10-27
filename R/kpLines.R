


kpLines <- function(karyoplot, data=NULL, chr=NULL, x=NULL, y=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, ...) {
  pp <- prepareParameters2("kpLines", karyoplot=karyoplot, data=data, chr=chr, x=x, y=y, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  ccf <- karyoplot$coord.change.function
  
  ss <- sapply(karyoplot$chromosomes, function(current.chr) {
    in.chr <- which(pp$chr==current.chr)
    xplot <- ccf(chr=pp$chr[in.chr], x=pp$x[in.chr], data.panel=data.panel)$x
    yplot <- ccf(chr=pp$chr[in.chr], y=pp$y[in.chr], data.panel=data.panel)$y
    lines(x=xplot, y=yplot, ...)      
  })  
}
