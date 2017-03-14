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
#' A valid \code{plot.params} object with the default values for the plotting parameters and 
#' ready to be used in the \code{plotKaryotype}
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
  valid.plot.types <- c(1:5) #c(1:4)
  
  if(!plot.type %in% valid.plot.types) {
    stop(paste0("plot.type is not valid. Select a valid value: ", 
                paste0(valid.plot.types, collapse=", ")))
  }
  
  if(plot.type == 1) { #Horizontal. Data above the ideogram
    plot.params <- list(leftmargin=0.1, rightmargin=0.05, topmargin=100, bottommargin=100,
                        ideogramheight=50, ideogramlateralmargin=0,
                        data1height=200, data1inmargin=20, data1outmargin=20,
                        data1min=0, data1max=1
    )
  } 
  if(plot.type == 2) { #Horizontal. Data above and below the ideogram
    plot.params <- list(leftmargin=0.1, rightmargin=0.05, topmargin=100, bottommargin=100,
                        ideogramheight=50, ideogramlateralmargin=0,
                        data1height=200, data1inmargin=20, data1outmargin=20,
                        data1min=0, data1max=1,
                        data2height=200, data2inmargin=20, data2outmargin=20,
                        data2min=0, data2max=1                        
                        )
  }
  if(plot.type == 3) { #Horizontal. All ideograms in a single line with 2 data panels
    plot.params <- list(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                        ideogramheight=10, ideogramlateralmargin=0,
                        data1height=200, data1inmargin=10, data1outmargin=0,
                        data1min=0, data1max=1,
                        data2height=200, data2inmargin=10, data2outmargin=0,
                        data2min=0, data2max=1                        
    )
  }
  if(plot.type == 4) { #Horizontal. All ideograms in a single line 1 panel above
    plot.params <- list(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                        ideogramheight=10, ideogramlateralmargin=0,
                        data1height=200, data1inmargin=10, data1outmargin=0,
                        data1min=0, data1max=1,
                        data2height=0, data2inmargin=0, data2outmargin=0, #make the panel invisible
                        data2min=0, data2max=1                        
    )
  }
  if(plot.type == 5) { #Horizontal. All ideograms in a single line 1 panel below
    plot.params <- list(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                        ideogramheight=10, ideogramlateralmargin=0,
                        data1height=0, data1inmargin=0, data1outmargin=0, #make this panel invisible
                        data1min=0, data1max=1,
                        data2height=200, data2inmargin=10, data2outmargin=0,
                        data2min=0, data2max=1                        
    )
  }
  
  
  return(plot.params)
}