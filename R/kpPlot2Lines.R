

#Given a color, returns a lighter one
#TODO: Is there a better way to do that?
lighter <- function(col, amount=150) {
  col <- (col2rgb(col)+amount)/255
  col[col[,1]>1,1] <- 1
  return(rgb(t(col)))  
}


#Helper Function
#From: http://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r
# segment-segment intersection code
# http://paulbourke.net/geometry/pointlineplane/
ssi <- function(x1, x2, x3, x4, y1, y2, y3, y4){
  
  denom <- ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1))
  denom[abs(denom) < 1e-10] <- NA # parallel lines
  
  ua <- ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom
  ub <- ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom
  
  x <- x1 + ua * (x2 - x1)
  y <- y1 + ua * (y2 - y1)
  inside <- (ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1)
  data.frame(x = ifelse(inside, x, NA), 
             y = ifelse(inside, y, NA))
  
}


kpPlot2Lines <- function(karyoplot, adata=NULL, achr=NULL, ax=NULL, ay=NULL, bdata=NULL, bchr=NULL, bx=NULL, by=NULL, ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, acol="red", bcol="blue", afill=NULL, bfill=NULL, ...) {

  app <- prepareParameters2("kpLines", karyoplot=karyoplot, data=adata, chr=achr, x=ax, y=ay, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  bpp <- prepareParameters2("kpLines", karyoplot=karyoplot, data=bdata, chr=bchr, x=bx, y=by, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
  
  ccf <- karyoplot$coord.change.function
  
  #Define the colors 
  if(is.null(afill)) afill <- lighter(acol)
  if(is.null(bfill)) bfill <- lighter(bcol)
  
  
  #Process separately for every chromosome
  
  for(chr in karyoplot$chromosomes) {
    #Filter A data
      a.in.chr <- app$chr == chr
      chr.achr <- app$chr[a.in.chr]
      chr.ax <- app$x[a.in.chr]
      chr.ay <- app$y[a.in.chr]
      a.len <- length(chr.achr)
    #Filter B data  
      b.in.chr <- bpp$chr == chr
      chr.bchr <- bpp$chr[b.in.chr]
      chr.bx <- bpp$x[b.in.chr]
      chr.by <- bpp$y[b.in.chr]
      b.len <- length(chr.bchr)
      
    #Find the intersection points between the 2 lines
    intersections <- ssi(x1=chr.ax[-a.len], x2=chr.ax[-1], y1=chr.ay[-a.len], y2=chr.ay[-1], #<- segments defined by line A
                         x3=chr.bx[-b.len], x4=chr.bx[-1], y3=chr.by[-b.len], y4=chr.by[-1]) # <- segments defined by line B
    intersections <- intersections[!is.na(intersections$x),]
  
    if(nrow(intersections)==0) { #If there are no intersections, simply plot a single polygon filing the space between the lines
        xpol <- c(chr.ax, rev(chr.bx))
        ypol <- c(chr.ay, rev(chr.by))
        #decide the color
        if(mean(chr.ay)>mean(chr.by)) {
          fillcol <- afill
        } else {
          fillcol <- bfill
        }
        #and plot the polygon
        kpPolygon(karyoplot=karyoplot, chr=chr, x=xpol, y=ypol, col=fillcol, border=NA, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)  
    } else {
      #If there are intersections, plot a polygon for every region they define
      #First, add a finall intersection point to plot the last polygon
      intersections <- rbind(intersections, c(x=mean(c(chr.ax[a.len], chr.bx[b.len])), y=mean(c(chr.ay[a.len], chr.by[b.len]))))
      #kpPoints(karyoplot, chr=chr, x=intersections$x, y=intersections$y)
          
      #These intersection points break the lines in regions where A is above B and regions where B is above A
      #Define these regions and plot the polygons defined by the 2 lines in the correct color
      old.ix <- NULL
      old.iy <- NULL
      for(n in c(1:nrow(intersections))) {
        ix <- intersections$x[n]
        iy <- intersections$y[n]
        if(is.null(old.ix)) { #If it's the first region
          a.in.region <- chr.ax < ix 
          b.in.region <- chr.bx < ix 
        } else { 
          a.in.region <- chr.ax <= ix & chr.ax >= old.ix
          b.in.region <- chr.bx <= ix & chr.ax >= old.ix
        }
        if(length(a.in.region)>0) { 
          #Prepare x coordinates
            xpol <- numeric()
            if(!is.null(old.ix)) xpol <- old.ix
            xpol <- c(xpol, chr.ax[a.in.region], ix, rev(chr.bx[b.in.region]))
          #Prepare y coordinates
            ypol <- numeric()
            if(!is.null(old.iy)) ypol <- old.iy
            ypol <- c(ypol, chr.ay[a.in.region], iy, rev(chr.by[b.in.region]))
          #decide the color
            if(mean(chr.ay[a.in.region])>mean(chr.by[b.in.region])) {
              fillcol <- afill
            } else {
              fillcol <- bfill
            }
          #and plot the polygon
            kpPolygon(karyoplot=karyoplot, chr=chr, x=xpol, y=ypol, col=fillcol, border=NA, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
        }
        old.ix <- ix
        old.iy <- iy
      }
    } 
    #And finally, plot the lines
      #Plot the line A
      kpLines(karyoplot=karyoplot, chr=app$chr, x=app$x, y=app$y, col=acol, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
    #Plot the line B
      kpLines(karyoplot=karyoplot, chr=bpp$chr, x=bpp$x, y=bpp$y, col=bcol, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
}
    
}
