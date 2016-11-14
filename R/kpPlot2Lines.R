# 
# 
# #Given a color, returns a lighter one
# #TODO: Is there a better way to do that?
# lighter <- function(col, amount=150) {
#   col <- (col2grDevices::rgbcol)+amount)/255
#   col[col[,1]>1,1] <- 1
#   return(grDevices::rgbt(col)))  
# }
# 
# #Helper Function
# #From: http://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r
# # segment-segment intersection code
# # http://paulbourke.net/geometry/pointlineplane/
# ssi <- function(x1, x2, x3, x4, y1, y2, y3, y4){
#   
#   denom <- ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1))
#   denom[abs(denom) < 1e-10] <- NA # parallel lines
#   
#   ua <- ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom
#   ub <- ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom
#   
#   x <- x1 + ua * (x2 - x1)
#   y <- y1 + ua * (y2 - y1)
#   inside <- (ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1)
#   data.frame(x = ifelse(inside, x, NA), 
#              y = ifelse(inside, y, NA))
#   
# }
# 
# 
# kpPlot2Lines <- function(karyoplot, adata=NULL, achr=NULL, ax=NULL, ay=NULL, 
#                          bdata=NULL, bchr=NULL, bx=NULL, by=NULL, 
#                          ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL, 
#                          acol="red", bcol="blue", afill=NULL, bfill=NULL, ...) {
# 
#   app <- prepareParameters2("kpPlot2Lines", karyoplot=karyoplot, data=adata, chr=achr, x=ax, y=ay, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel) , ...)
#   bpp <- prepareParameters2("kpPlot2Lines", karyoplot=karyoplot, data=bdata, chr=bchr, x=bx, y=by, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel) , ...)
#   
#   ccf <- karyoplot$coord.change.function
#   
#   #Define the colors 
#   if(is.null(afill)) afill <- lighter(acol, amount=200)
#   if(is.null(bfill)) bfill <- lighter(bcol, amount=200)
#   
# 
#   #Process separately for every chromosome
#   
#   aline <- data.frame(chr=app$chr, x=app$x, y=app$y, stringsAsFactors=FALSE)
#   bline <- data.frame(chr=bpp$chr, x=bpp$x, y=bpp$y, stringsAsFactors=FALSE)
# 
#   
#   for(chr in karyoplot$chromosomes) {
#     #Filter A data
#       a.in.chr <- app$chr == chr
#       chr.aline <- aline[a.in.chr,]
#       a.len <- nrow(chr.aline)
#     #Filter B data  
#       b.in.chr <- bpp$chr == chr
#       chr.bline <- bline[b.in.chr,]
#       b.len <- nrow(chr.bline)
#       
#     #Find the intersection points between the 2 lines
#     intersections <- ssi(x1=chr.aline$x[-a.len], x2=chr.aline$x[-1], y1=chr.aline$y[-a.len], y2=chr.aline$y[-1], #<- segments defined by line A
#                          x3=chr.bline$x[-b.len], x4=chr.bline$x[-1], y3=chr.bline$y[-b.len], y4=chr.bline$y[-1]) # <- segments defined by line B
#     intersections <- intersections[!is.na(intersections$x),]
#     intersections <- intersections[!duplicated(intersections$x),]
#       
#     if(nrow(intersections)==0) { #If there are no intersections, simply plot a single polygon filing the space between the lines
#         xpol <- c(chr.ax, rev(chr.bx))
#         ypol <- c(chr.ay, rev(chr.by))
#         #decide the color
#         if(mean(chr.ay)>mean(chr.by)) {
#           fillcol <- afill
#         } else {
#           fillcol <- bfill
#         }
#         #and plot the polygon
#         kpPolygon(karyoplot=karyoplot, chr=chr, x=xpol, y=ypol, col=fillcol, border=NA, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)  
#     } else {
#       #If there are intersections, plot a polygon for every region they define
#       na.a <- which(is.na(r$a$x) | is.nan(r$a$x) | is.na(r$a$y) | is.nan(r$a$y))
#       na.b <- which(is.na(r$b$x) | is.nan(r$b$x) | is.na(r$b$y) | is.nan(r$b$y))
#       cuts <- c(intersections$x, r$a$x[na.a], r$b$x[na.b])
#       cuts <- cuts[!is.na(cuts) & !is.nan(cuts)]
#       cuts <- sort(cuts)
#       cuts <- cuts[!duplicated(cuts)]
#       
#       #WARNING: We are ignoring the x=NA and x=NaN
#       
#       #break the lines using the cuts
#       ca <- chr.aline
#       cb <- chr.bline
#       
#       ca <- ca[-na.a,]
#       cb <- cb[-nb.b,]
#       
#       regions <- list()
#       for(x.cut in cuts) {
#         print(as.character(x.cut))
#         
#         #a line
#         a.before <- which(ca$x<=x.cut)
#         if(length(a.before)>0) {
#           reg.aline <- ca[a.before,]
#           print(paste0("ca: ", nrow(ca), "  a.before=", length(a.before), "   new ca:", nrow(ca[-a.before,])))
#           ca <- ca[-a.before,]
#         }
#         b.before <- which(cb$x<=x.cut)
#         if(length(b.before)>0) {
#           reg.bline <- ca[b.before,]
#           cb <- cb[-b.before,]
#         }
#         if(length(a.before)>0 | length(b.before)>0) {
#           regions[[as.character(x.cut)]] <- list(a=reg.aline, b=reg.bline)  
#         }
#       }
#       
#       
#       regions
#       
#       chr.aline$x
#       
#       
#       
#       
#       
#       #split the lines in contiguos regions using a) the intersection points b)the NA and NaN values
#       ca <- chr.aline
#       cb <- chr.bline
#       last.int <- NULL
#       regions <- list()
#       for(n in c(1:nrow(intersections))) {
#         ix <- intersections$x[n]
#         #a line
#           a.before <- which(ca$x<=ix)
#           if(length(a.before)>0) {
#             a.cut <- a.before[length(a.before)]
#             reg.aline <- ca[c(1:a.cut),]
#             ca <- ca[-c(1:a.cut),]
#           }
#         #a line
#           b.before <- which(cb$x<ix)
#           if(length(b.before)>0) {
#             b.cut <- b.before[length(b.before)]
#             reg.bline <- cb[c(1:b.cut),]
#             cb <- cb[-c(1:b.cut),]  
#           }
#         if(length(a.before)>0 | length(b.before)>0) {
#           regions[[as.character(n)]] <- list(a=reg.aline, b=reg.bline, ibefore=last.int, iafter=intersections[n,])  
#         }
#         last.int <- intersections[n,]
#       }
#       if(nrow(ca)>0 | nrow(cb)>0) {
#         regions[[as.character(n+1)]] <- list(a=ca, b=cb, ibefore=last.int, iafter=NULL)  
#       }
#       
#       #after that, split at NA and NaN to get the real regions to plot the polygons
#       final.regions <- list()
#       for(n in names(regions)) {
#         r <- regions[[n]]
#         na.a <- which(is.na(r$a$x) | is.nan(r$a$x) | is.na(r$a$y) | is.nan(r$a$y))
#         na.b <- which(is.na(r$b$x) | is.nan(r$b$x) | is.na(r$b$y) | is.nan(r$b$y))
#         if(length(na.a)>0 | length(na.b) > 0) { #If there are NAs or NaNs, split the region as needed
#           if(1 %in% na.a) { #if the first value of any of the lines is NA, forget the intersection point
#             na.a <- na.a[-1]
#             r$ibefore <- NULL
#           }
#           if(1 %in% na.b) {
#             na.b <- na.b[-1]
#             r$iafter <- NULL
#           }
#           
#           if(length(na.a)>0 | length(na.b) > 0) { #If there are still NAs
#             x.cuts <- sort(c(r$a$x[na.a-1], r$b$x[na.b-1]))
#             ra <- r$a[-na.a,]
#             rb <- r$b[-na.b,]        
#             
#             for(i in c(1:length(x.cuts))) {
#               x <- x.cuts[i]
#               if(is.na(x) | is.nan(x)) next; #do nothing with it
#               a.before <- which(ra$x <= x)
#               if(length(a.before)>0) {
#                 a <- ra[a.before,]
#                 ra <- ra[-a.before]
#               }
#               b.before <- which(rb$x <= x)
#               if(length(a.before)>0) {
#                 a <- ra[a.before,]
#                 ra <- ra[-a.before]
#               }
#               if(length(a)>1 & length(b)>1) { #If we cannot draw at least a 4 sided polygon, ignore the region and draw nothing (only lines)
#                                                 
#               } #else, do nothing. 
#             }  
#             
#             
#           } else {
#             final.regions[[n]] <- r
#           }
#         } else { #If there are no NAs or NaNs in the region, keep ot as is
#           final.regions[[n]] <- r
#         }
#         
#       }
#       
#       
#       n <- "8"
#       
#       #First, add a finall intersection point to plot the last polygon
#       intersections <- rbind(intersections, c(x=mean(c(chr.ax[a.len], chr.bx[b.len])), y=mean(c(chr.ay[a.len], chr.by[b.len]))))
#       #kpPoints(karyoplot, chr=chr, x=intersections$x, y=intersections$y)
#           
#       #These intersection points break the lines in regions where A is above B and regions where B is above A
#       #Define these regions and plot the polygons defined by the 2 lines in the correct color
#       old.ix <- NULL
#       old.iy <- NULL
#       for(n in c(1:nrow(intersections))) {
#         ix <- intersections$x[n]
#         iy <- intersections$y[n]
#         if(is.null(old.ix)) { #If it's the first region
#           a.in.region <- chr.ax < ix 
#           b.in.region <- chr.bx < ix 
#         } else { 
#           a.in.region <- chr.ax <= ix & chr.ax >= old.ix
#           b.in.region <- chr.bx <= ix & chr.ax >= old.ix
#         }
#         if(length(a.in.region)>0) { 
#           #Prepare x coordinates
#             xpol <- numeric()
#             if(!is.null(old.ix)) xpol <- old.ix
#             xpol <- c(xpol, chr.ax[a.in.region], ix, rev(chr.bx[b.in.region]))
#           #Prepare y coordinates
#             ypol <- numeric()
#             if(!is.null(old.iy)) ypol <- old.iy
#             ypol <- c(ypol, chr.ay[a.in.region], iy, rev(chr.by[b.in.region]))
#           #decide the color
#             if(mean(chr.ay[a.in.region])>mean(chr.by[b.in.region])) {
#               fillcol <- afill
#             } else {
#               fillcol <- bfill
#             }
#           #and plot the polygon
#             kpPolygon(karyoplot=karyoplot, chr=chr, x=xpol, y=ypol, col=fillcol, border=NA, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel)
#         }
#         old.ix <- ix
#         old.iy <- iy
#       }
#     } 
#     #And finally, plot the lines
#       #Plot the line A
#       kpLines(karyoplot=karyoplot, chr=app$chr, x=app$x, y=app$y, col=acol, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
#     #Plot the line B
#       kpLines(karyoplot=karyoplot, chr=bpp$chr, x=bpp$x, y=bpp$y, col=bcol, ymin=ymin, ymax=ymax, r0=r0, r1=r1, data.panel=data.panel, ...)
#   }
#     
# }
