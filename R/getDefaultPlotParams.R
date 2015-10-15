


getDefaultPlotParams <- function(plot.type) {
  valid.plot.types <- c(1:4)
  
  if(!plot.type %in% valid.plot.types) {
    stop(paste0("plot.type is not valid. Select a valid value: ", paste0(valid.plot.types, collapse=", ")))
  }
  
  if(plot.type == 1) { #Horizontal. Data above the ideogram
    plot.params <- list(xleftmargin=0.1, xrightmargin=0.05, ytopmargin=100, ybottommargin=100,
                      yabovemargin=10, ybelowmargin=10, ydataheight=200, ideogramheight=50,
                      dataymin=0, dataymax=1)
  } 
  if(plot.type == 2) { #Horizontal. Data above and below the ideogram
    plot.params <- list(xleftmargin=0.1, xrightmargin=0.05, ytopmargin=100, ybottommargin=100,
                        yabovemargin=10, ybelowmargin=10, ydataheight=200, ideogramheight=50,
                        dataymin=0, dataymax=1)
  }
  if(plot.type == 3) { #Horizontal. Data below the ideogram
    stop("Plot type 3 still unimplemented")
  }
  if(plot.type == 4) { #Horizontal Linear. Data above the ideogram.
    stop("Plot type 4 still unimplemented")
  }
  
  return(plot.params)
}