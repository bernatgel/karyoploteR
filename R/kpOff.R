
#' @export kpOff

#Finalize the KaryPlot (basically, recover the prior graphical parameters)

kpOff <- function(karyoplot) {
  if(!is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  par(karyoplot$plot$old.par)
    
}