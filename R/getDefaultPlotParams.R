#' getDefaultParameters
#' 
#' @description 
#' 
#' Returns the default parameters for the given plot.type
#' 
#' @details 
#'  
#'  Given a plot.type, this function returns a list suitable as a valid \code{plot.params} object.
#'  The user can then proceed to change the parameter values as needed and supply the modified 
#'  list to the plotKaryotype function.#'
#'  
#' @usage getDefaultPlotParams(plot.type)
#'  
#' @param plot.type   (integer) the required plot type. can be any valid plot type (see \code{\link{plotKaryotype}})
#'
#' @return 
#' 
#' A valid \code{plot.params} object with the default values for the plotting parameters and ready to be used in the \code{plotKaryotype}
#'
#' @seealso \code{\link{plotKaryotype}}
#' 
#' @examples
#' 
#' pp <- getDefaultPlotParams(plot.type=2)
#' pp
#'
#' #Change the ideogramheight param to create thicker ideograms 
#' pp$ideogramheight <- 150
#' 
#' plotKaryotype(genome="hg19", plot.type=2, plot.params=pp) 
#' 
#'  
#' @export getDefaultPlotParams
#' 


getDefaultPlotParams <- function(plot.type) {
  valid.plot.types <- c(1:2) #c(1:4)
  
  if(!plot.type %in% valid.plot.types) {
    stop(paste0("plot.type is not valid. Select a valid value: ", paste0(valid.plot.types, collapse=", ")))
  }
  
  if(plot.type == 1) { #Horizontal. Data above the ideogram
    plot.params <- list(leftmargin=0.1, rightmargin=0.05, topmargin=100, bottommargin=100,
                        ideogramheight=50, 
                        data1height=200, data1inmargin=20, data1outmargin=20, data1min=0, data1max=1
    )
  } 
  if(plot.type == 2) { #Horizontal. Data above and below the ideogram
    plot.params <- list(leftmargin=0.1, rightmargin=0.05, topmargin=100, bottommargin=100,
                        ideogramheight=50, 
                        data1height=200, data1inmargin=20, data1outmargin=20, data1min=0, data1max=1,
                        data2height=200, data2inmargin=20, data2outmargin=20, data2min=0, data2max=1                        
                        )
  }
  if(plot.type == 3) { #Horizontal. Data below the ideogram
    stop("Plot type 3 still unimplemented")
  }
  if(plot.type == 4) { #Horizontal Linear. Data above the ideogram.
    plot.params <- list(leftmargin=0.1, rightmargin=0.05, topmargin=100, bottommargin=100,
                        ideogramheight=50, 
                        data1height=200, data1inmargin=20, data1outmargin=20, data1min=0, data1max=1,
                        data2height=200, data2inmargin=20, data2outmargin=20, data2min=0, data2max=1                        
    )
  }
  
  return(plot.params)
}