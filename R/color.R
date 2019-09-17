############  Colors  ###############

############ Constants ##############

#This structure contains the different color schemas used in karyoploteR
.karyoploter.colors <- list(
  cytobands=list(
    schemas=list(
      circos=c(gneg="#FFFFFF",
               gpos25="#C8C8C8",
               gpos33="#D2D2D2",
               gpos50="#C8C8C8",
               gpos66="#A0A0A0",
               gpos75="#828282",
               gpos100="#000000",
               gpos="#000000",
               stalk="#647FA4", #repetitive areas
               acen="#D92F27", #centromeres
               gvar="#DCDCDC",
               border="black"),
      only.centromeres=c(gneg="#C8C8C8",
                         gpos25="#C8C8C8",
                         gpos33="#C8C8C8",
                         gpos50="#C8C8C8",
                         gpos66="#C8C8C8",
                         gpos75="#C8C8C8",
                         gpos100="#C8C8C8",
                         gpos="#C8C8C8",
                         stalk="#C8C8C8", #repetitive areas
                         acen="#D92F27", #centromeres
                         gvar="#C8C8C8",
                         border=NA),
      biovizbase=c(gneg = "grey100", stalk = "brown3", acen = "brown4", gpos = "grey0", 
                   gvar = "grey0", gpos1 = "#FFFFFF", gpos2 = "#FCFCFC", gpos3 = "#F9F9F9", 
                   gpos4 = "#F7F7F7", gpos5 = "#F4F4F4", gpos6 = "#F2F2F2", gpos7 = "#EFEFEF", 
                   gpos8 = "#ECECEC", gpos9 = "#EAEAEA", gpos10 = "#E7E7E7", gpos11 = "#E5E5E5", 
                   gpos12 = "#E2E2E2", gpos13 = "#E0E0E0", gpos14 = "#DDDDDD", gpos15 = "#DADADA", 
                   gpos16 = "#D8D8D8", gpos17 = "#D5D5D5", gpos18 = "#D3D3D3", gpos19 = "#D0D0D0", 
                   gpos20 = "#CECECE", gpos21 = "#CBCBCB", gpos22 = "#C8C8C8", gpos23 = "#C6C6C6", 
                   gpos24 = "#C3C3C3", gpos25 = "#C1C1C1", gpos26 = "#BEBEBE", gpos27 = "#BCBCBC", 
                   gpos28 = "#B9B9B9", gpos29 = "#B6B6B6", gpos30 = "#B4B4B4", gpos31 = "#B1B1B1", 
                   gpos32 = "#AFAFAF", gpos33 = "#ACACAC", gpos34 = "#AAAAAA", gpos35 = "#A7A7A7", 
                   gpos36 = "#A4A4A4", gpos37 = "#A2A2A2", gpos38 = "#9F9F9F", gpos39 = "#9D9D9D", 
                   gpos40 = "#9A9A9A", gpos41 = "#979797", gpos42 = "#959595", gpos43 = "#929292", 
                   gpos44 = "#909090", gpos45 = "#8D8D8D", gpos46 = "#8B8B8B", gpos47 = "#888888", 
                   gpos48 = "#858585", gpos49 = "#838383", gpos50 = "#808080", gpos51 = "#7E7E7E", 
                   gpos52 = "#7B7B7B", gpos53 = "#797979", gpos54 = "#767676", gpos55 = "#737373", 
                   gpos56 = "#717171", gpos57 = "#6E6E6E", gpos58 = "#6C6C6C", gpos59 = "#696969", 
                   gpos60 = "#676767", gpos61 = "#646464", gpos62 = "#616161", gpos63 = "#5F5F5F", 
                   gpos64 = "#5C5C5C", gpos65 = "#5A5A5A", gpos66 = "#575757", gpos67 = "#545454", 
                   gpos68 = "#525252", gpos69 = "#4F4F4F", gpos70 = "#4D4D4D", gpos71 = "#4A4A4A", 
                   gpos72 = "#484848", gpos73 = "#454545", gpos74 = "#424242", gpos75 = "#404040", 
                   gpos76 = "#3D3D3D", gpos77 = "#3B3B3B", gpos78 = "#383838", gpos79 = "#363636", 
                   gpos80 = "#333333", gpos81 = "#303030", gpos82 = "#2E2E2E", gpos83 = "#2B2B2B", 
                   gpos84 = "#292929", gpos85 = "#262626", gpos86 = "#242424", gpos87 = "#212121", 
                   gpos88 = "#1E1E1E", gpos89 = "#1C1C1C", gpos90 = "#191919", gpos91 = "#171717", 
                   gpos92 = "#141414", gpos93 = "#121212", gpos94 = "#0F0F0F", gpos95 = "#0C0C0C", 
                   gpos96 = "#0A0A0A", gpos97 = "#070707", gpos98 = "#050505", gpos99 = "#020202", 
                   gpos100 = "#000000", border = "black")
    )
  ),
  variants=list(
    schemas=list(
      cell21breast=c("C>A"="#4c64ae",
        "C>G"="#000000",
        "C>T"="#e40611",
        "T>A"="#bf4a96",
        "T>C"="#fbe800",
        "T>G"="#6eb529",
        "other"="#888888")
    )
  ),
  horizon=list(
    schemas=list(
      redblue6=c("#b51010", "#e73131", "#f7ad9c", "#ffffff", "#a5c6de", "#4284c6", "#005aad"),
      bluepurple10=c("#CAE5FF", "#AAD6FF", "#89C6FF", "#AAD6FF", "#CAE5FF", "#AFD2E9", "#9D96B8", "#9a7197", "#886176", "#7C5869"),
      bluegold3=c("#4284c6", "white", "gold")
    )
  )
)


#' getColorSchemas
#' 
#' @description 
#' Return a structure with the color schemas included in karyoploteR
#' 
#' @usage getColorSchemas()
#' 
#' @return
#' A list with the color schemas included in karyoploteR for cytobands,
#' variants, horizons...
#' 
#' @examples
#'  
#' getColorSchemas()
#'  
#' @export is.color
getColorSchemas <- function() {
  return(.karyoploter.colors)
}


#' plotPalette
#' 
#' @description 
#' Create a plot of the palette.
#' 
#' @details 
#' Creates a simple plot with a rectangle for every color of every palette.
#' cols must be either a vector of colors (in any format accepted by 
#' karyoploteR::"is.color") or a list of such vectors.
#' The names of the list elements will be treated as the palette names, if
#' the list has no names, palettes will be called "Pallete1", "Palette2", ...
#' 
#' @usage plotPalettes(cols, add.color.name=TRUE, border=NA, palette.names.col="black", palette.names.cex=1, palette.names.srt=0, color.names.col="black", color.names.cex=1, color.names.srt=0, ...)
#' 
#' @param cols (color vector or list of color vectors) The colors to plot
#' @param add.color.name (logical) Wether to add or not the names of the colors, their definition.
#' @param palette.names.col (color) The color of the palette names (defaults to "black")
#' @param palette.names.cex (numeric) The cex value (size) for the palette names (defaults to 1)
#' @param palette.names.srt (numeric) The srt value (rotation) for the palette names (defaults to 0)
#' @param color.names.col (color) The color of the color names (defaults to "black")
#' @param color.names.cex (numeric) The cex value (size) for the color names  (defaults to 1)
#' @param color.names.srt (numeric) The srt value (rotation) for the color names (defaults to 0)
#' @param border (color) The color of the border of the palette rectangles. If NA, no border. (defaults to NA)
#' @param ... Any additional plotting parameters
#' 
#' @return
#' nothing
#' 
#' 
#' @examples
#'  
#' plotPalettes(c("red", "blue", "yellow", "green"))
#' palettes <- list("P1"=c("red", "#000000", lighter("gold")), 
#'                  "P2"=c("orchid", "yellow"))
#' plotPalettes(palettes, color.names.col=c("blue", "green", "red"), border="black", color.names.srt=45)
#' 
#' @export plotPalettes
#' 
plotPalettes <- function(cols, add.color.name=TRUE, border=NA, palette.names.col="black", palette.names.cex=1, palette.names.srt=0, color.names.col="black", color.names.cex=1, color.names.srt=0, ...) {
  if(!is.list(cols)) {
    cols <- list(cols)
  }
  if(any(!unlist(lapply(cols, is.color)))) stop("All elements in cols must be valid colors.")
  
  #Set the palette names if not available
  if(is.null(names(cols))) names(cols) <- paste0("Palette", seq_along(cols))
  
  #create an empty plot
  graphics::plot(x=0, type="n", 
       xlim=c(0, 10*max(unlist(lapply(cols, length)))), ylim=c(0, 10*length(cols)),
       ylab="", xlab="", axes=FALSE, xaxs="i", yaxs="i")
  
  #Add the names of the palettes as y labels
  graphics::axis(2, at = 4+(seq_along(cols)-1)*10, labels=rev(names(cols)), las=2, 
                 lwd = 0, cex=palette.names.cex, col=palette.names.col, srt=palette.names.srt)   
  
  #Plot the palettes as rectangles
  for(npal in seq_along(cols)) {
    pp <- cols[[length(cols) - npal + 1]]
    graphics::rect(xleft=10*(seq_along(pp)-1), xright=8+10*(seq_along(pp)-1), 
         ybottom = (npal-1)*10, ytop = 8+(npal-1)*10, col=pp, border = border, ...)
    if(add.color.name) {
      graphics::text(x=4+10*(seq_along(pp)-1), y=4+(npal-1)*10, 
                     labels=as.character(pp), col=color.names.col,
                     cex=color.names.cex, srt=color.names.srt)
    }
  }
}

#' lighter
#' 
#' @description 
#' Given a color, return a lighter one
#' 
#' @details 
#' Very simple utility function to create lighter colors. Given a color, it
#' transforms it to rgb space, adds a set amount to all chanels and transforms
#' it back to a color.
#' 
#' @usage lighter(col, amount=150)
#' 
#' @param col (color) The original color. Might be specified as a color name or a "#RRGGBB(AA)" hex color definition.
#' @param amount (integer, [0-255]) The fixed amount to add to each RGB channel (Defaults to 150).
#' 
#' @return
#' A lighter color
#' 
#' @seealso \code{\link{darker}}
#' 
#' @examples
#'  
#' lighter("red")
#' lighter("#333333")
#' lighter(c("red", 3, "#FF00FF"))
#'  
#' @export lighter
#' 

lighter <- function(col, amount=150) {
  #Colors must be specified by name or #RRGGBB(AA)
  if(!all(is.color(col))) stop("All elements in col must be valid colors. Use is.col(col) to check it.")
  if(!methods::is(amount, "numeric") || length(amount)!=1) stop("amount must be a single number")
  if(amount>255 || amount<0) stop("amount must be a number between 0 and 255")
  .lighter <- function(col, amount) {
    if(is.na(col)) return(NA)
    new.col <- ((grDevices::col2rgb(col))+amount)/255
    new.col[new.col[,1]>1,1] <- 1
    return(grDevices::rgb(t(new.col)))  
  }
  
  return(unlist(lapply(col, .lighter, amount)))
}

#' darker
#' 
#' @description 
#' Given a color, return a darker one
#' 
#' @details 
#' Very simple utility function to create darker colors. Given a color, it
#' transforms it to rgb space, adds a set amount to all chanels and transforms
#' it back to a color.
#' 
#' @usage darker(col, amount=150)
#' 
#' @param col (color) The original color. Might be specified as a color name or a "#RRGGBB(AA)" hex color definition.
#' @param amount (integer, [0-255]) The fixed amount to subtract to each RGB channel (Defaults to 150).
#' 
#' @return
#' A darker color
#' 
#' @seealso \code{\link{lighter}}
#' 
#' @examples
#'  
#' darker("red")
#' darker("#333333")
#' darker(c("red", 3, "#FF00FF"))
#'  
#' @export darker
#'

#Given a color, returns a darker one
darker <- function(col, amount=150) {
  #Colors must be specified by name or #RRGGBB(AA)
  if(!all(is.color(col))) stop("All elements in col must be valid colors. Use is.col(col) to check it.")
  if(!methods::is(amount, "numeric") || length(amount)!=1) stop("amount must be a single number")
  if(amount>255 || amount<0) stop("amount must be a number between 0 and 255")
  
  .darker <- function(col, amount) {
    if(is.na(col)) return(NA)
    new.col <- ((grDevices::col2rgb(col))-amount)/255
    new.col[new.col[,1]<0, 1] <- 0
    return(grDevices::rgb(t(new.col)))
  }
  
  return(unlist(lapply(col, .darker, amount)))
}


#' transparent
#' 
#' @description 
#' Given a color, return a transparent one
#' 
#' @details 
#' Very simple utility function to create transparent colors. Given a color, it
#' transforms it to rgb space, adds a set amount to all chanels and transforms
#' it back to a color.
#' 
#' @usage transparent(col, amount=0.5)
#' 
#' @param col (color) The original color. Might be specified as a color name or a "#RRGGBB(AA)" hex color definition.
#' @param amount (number, [0-1]) The amount of transparency. 0 for completely visible, 1 for completely transparent. (Defaults to 0.5).
#' 
#' @return
#' A transparent color
#' 
#' @seealso \code{\link{lighter}}
#' 
#' @examples
#'  
#' transparent("red")
#' transparent("#333333")
#'  
#' @export transparent
#' @importFrom grDevices adjustcolor

#Given a color, returns a transparent one
transparent <- function(col, amount=0.5) {
  #Colors must be specified by name or #RRGGBB(AA)
  if(!methods::is(col, "character")) stop("Unknown color definition.")
  if(!methods::is(amount, "numeric") || length(amount)!=1) stop("amount must be a single number")
  if(amount>1 || amount<0) stop("amount must be a number between 0 and 1")
  
  transp.col <- grDevices::adjustcolor(col, alpha.f = (1-amount))
  transp.col[is.na(col)] <- NA
  return(transp.col)
}



#' colByChr
#' 
#' @description 
#' Given a set of data elements, return a color for each one based on their chromosome
#' 
#' @details 
#' Returns a color for each data element based on its chromosome. The returned colors might
#' com from one of the predefined color sets or passed in as a parameter.
#' 
#' If \code{colors} is the name of one of the available color sets, it the color set is used. 
#' If it's a named character vector with the chromosome as names, they will be assigned by name
#' and any missing chromosome will be \code{default.col}. If it's a non-named chraracter vector,
#' will be used in order and recycled if necessary.
#' 
#' Data might be either a GRanges object or a vector of chromosomes.
#' 
#' @usage colByChr(data, colors="2grays", all.chrs=NULL, default.col="black")
#' 
#' @param data Either a vector of characters or a GRanges object
#' @param colors The name of a color set ("2grays", "blackgreen", "rainbow"...) or a vector of colors. If the vector is named, names are expected to be the chromosome names. (defaults to "2grays")
#' @param all.chrs A vector with all possible chromosomes. If NULL, the list will be extracted from data (using seqlevels if available). (defaults to NULL)
#' @param default.col The default color to return when something is unavailable
#' 
#'
#' @return
#' A vector of colors
#'
#' @note 
#' Available color.sets:
#' "2grays"=c("#888888", "#444444"),
#' "2blues"=c("#6caeff", "#2b5d9b")
#' "blackgreen"=c("black", "green"),
#' "greengray"=c("#c6ffb7", "#888888"),
#' "brewer.set1"=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
#' "brewer.set2"=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
#' "brewer.set3"=c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
#' "brewer.pastel1"=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"),
#' "brewer.pastel2"=c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC"),
#' "rainbow"=rainbow(n=length(all.chrs))
#'   
#' @seealso \code{\link{kpPoints}}
#' 
#' @examples
#' 
#' chrs <- c("chr1", "chr2", "chr2", "chr1", "chr5")
#' points <- toGRanges(paste0("chr", c(1:22, "X", "Y")), rep(10e6, 24), rep(10e6, 24))
#' 
#' colByChr(chrs)
#' colByChr(points)
#' 
#' kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, ideogram.plotter=NULL)
#' kpAddChromosomeNames(kp, srt=45)
#' kpAddChromosomeSeparators(kp)
#' 
#' total.tracks <- 6
#' 
#' kpPoints(kp, points, col=colByChr(points), y=0.5, cex=1, r0=autotrack(1,total.tracks)$r0, r1=autotrack(1,total.tracks)$r1)
#' colors <- NULL
#' kpPoints(kp, points, y=0.5, col=colByChr(points, colors=colors), cex=1, r0=autotrack(2,total.tracks)$r0, r1=autotrack(2,total.tracks)$r1)
#' colors <- c("red", "blue")
#' kpPoints(kp, points, y=0.5, col=colByChr(points, colors=colors), cex=1, r0=autotrack(3,total.tracks)$r0, r1=autotrack(3,total.tracks)$r1)
#' colors <- c(chr1="red", chr7="blue")
#' kpPoints(kp, points, y=0.5, col=colByChr(points, colors=colors), cex=1, r0=autotrack(4,total.tracks)$r0, r1=autotrack(4,total.tracks)$r1)
#' kpPoints(kp, points, y=0.5, col=colByChr(points, colors=colors, default.col="green"), cex=1, r0=autotrack(5,total.tracks)$r0, r1=autotrack(5,total.tracks)$r1)
#' colors <- c("red", "yellow", 3, "orchid", "blue")
#' kpPoints(kp, points, y=0.5, col=colByChr(points, colors=colors), cex=1, r0=autotrack(6,total.tracks)$r0, r1=autotrack(6,total.tracks)$r1)
#' 
#' #Color sets
#' pp <- getDefaultPlotParams(plot.type=4)
#' pp$leftmargin <- 0.2
#' kp <- plotKaryotype(plot.type=4, labels.plotter=NULL, ideogram.plotter=NULL, plot.params=pp)
#' kpAddChromosomeNames(kp, srt=45)
#' kpAddChromosomeSeparators(kp)
#' 
#' color.sets <- c( "2grays", "2blues", "blackgreen", "greengray", "brewer.set1",
#'                    "brewer.set2", "brewer.set3", "brewer.pastel1", "brewer.pastel2", "rainbow" )
#' total.tracks <- length(color.sets)
#' for(i in seq_len(length(color.sets))) {
#'     kpPoints(kp, points, y=0.5, col=colByChr(points, colors=color.sets[i]), cex=1, r0=autotrack(i,total.tracks)$r0, r1=autotrack(i,total.tracks)$r1)
#'     kpAddLabels(kp, labels=color.sets[i], cex=0.7, r0=autotrack(i,total.tracks)$r0, r1=autotrack(i,total.tracks)$r1)
#' }
#' 
#' @export colByChr
#' @importFrom grDevices rainbow
#'

colByChr <- function(data, colors="2grays", all.chrs=NULL, default.col="black") {
  #Process data
  chrs <- NULL
  if(is.character(data)) {
    chrs <- data
  } else { #It's not a vector of chromosomes, try to extract them for other typical objects
    #GRanges
    chrs <- tryCatch(as.character(seqnames(data)), error=function(e) {return(NULL)})
    #TODO: VCF (we'll need to odepend on variant annotation)
    #if(is.null(chrs)) chrs <- tryCatch(as.character(seqnames(rowRanges(data))), error=function(e) {return(NULL)})
    #TODO: Summarized Experiment?
  }
  
  if(is.null(chrs)) {
    stop("Unknown data type passed to colByChr.")
  }
  
  #Try to get the chromosomes "in the right order" if possible
  if(is.null(all.chrs)) {
    all.chrs <- tryCatch(as.character(seqlevels(data)), error=function(e) {return(NULL)})
    if(is.null(all.chrs)) { #If it's not been possible, just get them in any order
      all.chrs <- unique(chrs)
    }
  }
  
  #Process color
  cols <- NULL
  if(!is.null(colors)) {
    color.sets <- list(
      "2grays"=c("#888888", "#444444"),
      "2blues"=c("#6caeff", "#2b5d9b"),
      "blackgreen"=c("black", "green"),
      "greengray"=c("#c6ffb7", "#888888"),
      "brewer.set1"=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
      "brewer.set2"=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
      "brewer.set3"=c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"),
      "brewer.pastel1"=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2"),
      "brewer.pastel2"=c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC"),
      "rainbow"=rainbow(n=length(all.chrs))
    )
    
    if(is.character(colors) && length(colors)==1L && colors %in% names(color.sets)) { #Name of a color set
      colors <- color.sets[[colors]]
    } 
    
    #Named color vector
    if(!is.null(names(colors))) {
      missing.chrs <- all.chrs[!(all.chrs %in% names(colors))] 
      cols <- c(colors, setNames(rep(default.col, length(missing.chrs)), missing.chrs))
    } else {
      cols <- setNames(rep(colors, length = length(all.chrs)), all.chrs)
    }   
  }
  if(is.null(cols)) cols <- setNames(rep(default.col, length(all.chrs)), all.chrs)
  
  #We now have chrs and cols, return thee colors per chromosome
  return(cols[chrs])
}




#' colByRegion
#' 
#' @description 
#' Given a set of data elements, return a color for each one based on whether 
#' they overlap a given set of regions. This might be useful, for example,  to
#' set a different color for data points overlapping a certain region of
#' interest.
#' 
#' @details 
#' Given a set of data elements, return a color for each one based on whether 
#' they overlap a given set of regions. The colors might be different for each
#' region and can be specified either in the regions object itself or in 
#' a separate \code{colors} parameter. If specified in \code{colors}, the values
#' will be recycled as needed. 
#' 
#' 
#' @usage colByRegion(data, regions, colors=NULL, default.col="black")
#' 
#' @param data Either a vector of characters or a GRanges object
#' @param regions (GRanges or equivalent) A set of regions where the color will be modified. Internally it will be converted into a Genomic Ranges object by \code{\link[regioneR]{toGRanges}} (from regioneR package) and so it can be either a GRanges, a data.frame, a character or any other value type accepted by that function. If \code{colors} is NULL (the default) and regions has additional columns in addition to chr, start and end, if any has a name in c("color", "colors", "col", "cols") it will be used. Otherwise the first additional column will be used.
#' @param colors (color) The colors to be used for each region. The content will be recycled if needed. If NULL, the colors are assumed to be available in the regions object. (defaults to NULL)
#' @param default.col The default color to return for data elements not overlapping the regions. (defaults to "black")
#' 
#'
#' @return
#' A vector of colors
#'
#' @seealso \link{kpPoints}, \link{colByChr}, \link[regioneR]{toGRanges}
#' 
#' @examples
#' 
#' data <- toGRanges("chr1", c(1e6*1:245), c(1e6*1:245)+10)
#' data$y <- rnorm(n = length(data), mean = 0.5, sd = 0.15)
#' 
#' regions <- toGRanges(c("chr1:10e6-20e6", "chr1:100e6-150e6"))
#' regions$col <- c("red", "blue")
#' 
#' kp <- plotKaryotype(chromosomes="chr1")
#' kpPoints(kp, data=data, r0=0, r1=0.2)
#' kpPoints(kp, data=data, r0=0.2, r1=0.4, col=colByRegion(data, regions = regions) )
#' kpText(kp, data=data, r0=0.4, r1=0.6, col=colByRegion(data, regions = regions), label="A", cex=0.5 )
#' kpBars(kp, data=data, y0=0, y1=data$y, r0=0.6, r1=0.8, border=colByRegion(data, regions = regions))
#' #It might not work wor objects where R expects a single color such as lines. Segments should be used instead
#' kpLines(kp, data=data, r0=0.8, r1=1, col=colByRegion(data, regions = regions) )
#' 
#' 
#' kp <- plotKaryotype(chromosomes="chr1")
#' kpPoints(kp, data=data, r0=0, r1=0.25)
#' kpPoints(kp, data=data, r0=0.25, r1=0.5, col=colByRegion(data, regions = regions, colors="green") )
#' kpText(kp, data=data, r0=0.5, r1=0.75, col=colByRegion(data, regions = regions, color=c("gray", "gold")), label="A", cex=0.5 )
#' kpBars(kp, data=data, y0=0, y1=data$y, r0=0.75, r1=1, border=colByRegion(data, regions = regions))
#' 
#' @export colByRegion
#'
#' @importFrom  GenomeInfoDb seqlevelsStyle

colByRegion <- function(data, regions, colors=NULL, default.col="black") {
  data <- toGRanges(data)
  
  regions <- toGRanges(regions)
  
  #Get the colors if needed
  if(is.null(colors)) {
    #Get the color from the regions GRanges
    if(length(mcols(regions))>0) {
      colors.column <- which(names(mcols(regions)) %in% c("color", "colors", "col", "cols"))[1]
      #if no column has a valid name, use the first one
      if(is.na(colors.column)) colors.column <- 1
      colors <- mcols(regions)[,colors.column]
    }
  } else {
    colors <- rep(colors, length.out=length(regions))
  }
  
  if(!all(is.color(colors))) { stop("A valid color specification is needed. Either as a column of regions or in the 'colors' parameter") }
  
  #Make regions and data comparable
  #WARNING: We are doing this here but nowhere else in the package. Should we?
  GenomeInfoDb::seqlevelsStyle(regions) <- GenomeInfoDb::seqlevelsStyle(data)
  
  #Assign everything the default color and then change it for the data
  #overlapping the regions
  cols <- rep(default.col, length(data))
  for(num.reg in seq_len(length(regions))) {
    cols[overlapsAny(data, regions[num.reg])] <- colors[num.reg]
  }
  
  return(cols)
}




#' is.color
#' 
#' @description 
#' Test if something is a valid color
#' 
#' @details 
#' This function tests if something is a valid color. Returns TRUE or FALSE. 
#' The function is vectorised. 
#' 
#' 
#' @usage is.color(x)
#' 
#' @param x The element to test
#' 
#' @return
#' TRUE is x is a valid color, FALSE otherwise
#' 
#' 
#' @examples
#'  
#' is.color("red")
#' is.color("#333333")
#' is.color(NA)
#' is.color(NULL)
#' is.color("not_a_color")
#' is.color(3)
#' 
#' is.color(c("not_a_color", "red", 3, "#FF0000"))
#'  
#' @export is.color
#' 
#' @importFrom grDevices col2rgb
#' 
#'
# Adapted from https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation

is.color <- function(x) {
  if(is.null(x)) return(FALSE)
  return(setNames(vapply(x, function(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)), 
             error = function(e) FALSE)
  },
  FUN.VALUE = TRUE), NULL))
}






#TODO: Document and export?

#Preprocess col and border to assign one based on the other if any is NULL
preprocessColors <- function(col=NULL, border=NULL, default.col="gray70", amount=100) {
  #If both NULL, return the default values
  if(is.null(col) && is.null(border)) {
    return(list(col=default.col, border=darker(default.col, amount = amount)))
  }
  
  if(is.null(border)) {
    if(all(is.na(col))) {
      border <- darker(default.col, amount = amount)
    } else {
      border <- darker(col, amount = amount)
    }
  }
  
  if(is.null(col)) {
    if(all(is.na(border))) {
      col <- default.col
    } else {
      col <- lighter(border, amount = amount)
    }
  }
  

  return(list(col=col, border=border))
}








############################# Function specific color seleectors ##############################
#' horizonColors
#' 
#' @description 
#' Returns the color structure needed by kpPlotHorizon
#' 
#' @details 
#' This function transforms an array of colors into a list of colors 
#' **internally** needed by kpPlotHorizon: a list with two elements, "neg" 
#' and "pos", each an array of colors of length num.parts. If col is 
#' a character of length one, it is interpreted as the name of a color scheme.
#' 
#' 
#' horizonColors(col, num.parts) 
#' 
#' @param col (array of colors) An array of colors
#' @param num.parts (positive integer) The number of colors to generate for pos and neg
#' 
#' @return
#' A list with 2 elements, pos and neg, each with num.parts colors
#' 
#' 
#' @examples
#'  
#' horizonColors("redblue6", 3)
#' horizonColors("redblue6", 6)
#' horizonColors("bluegold3", 2)
#' horizonColors(c("red", "blue"), 3)
#' horizonColors(c("red", "#FFFFFF00", "blue"), 3)
#'  
#' @export horizonColors
#' @importFrom grDevices colorRamp
horizonColors <- function(col, num.parts) {
  if(is.character(col) && length(col)==1) col <- .karyoploter.colors$horizon$schemas[[col]]
  ramp <- grDevices::colorRamp(col, alpha=TRUE)
  num.cols <- num.parts*2
  cols <- ramp(1/(num.cols)*c(0:num.cols))/255
  cols <- grDevices::rgb(cols, alpha = cols[,4])
  return(list(neg=rev(cols[1:num.parts]),
              pos=cols[(num.parts+2):(2*num.parts+1)]
  ))
}
