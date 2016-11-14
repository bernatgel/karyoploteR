# #INTERNAL
# 
# #Returns a table with the colors associated with each cytoband.
# # By default it returns the colors defined in Circos.
# 
# getCytobandColors <- function(color.table) {
#   
#   if(is.null(color.table)) {
#     color.table <- list(gneg="#FFFFFF",
#                      gpos25="#C8C8C8",
#                      gpos33="#D2D2D2",
#                      gpos50="#C8C8C8",
#                      gpos66="#A0A0A0",
#                      gpos75="#828282",
#                      gpos100="#000000",
#                      gpos="#000000",
#                      stalk="#647FA4", #repetitive areas
#                      acen="#D92F27", #centromeres
#                      gvar="#DCDCDC")
#   }
#   return(color.table)
# }
# 
# # Default UCSC color scheme for chromosome colors. 
# 
# chr1  = grDevices::rgb153/256,102/256,0/256)
# chr2  = 102,102,0
# chr3  = 153,153,30
# chr4  = 204,0,0
# chr5  = 255,0,0
# chr6  = 255,0,204
# chr7  = 255,204,204
# chr8  = 255,153,0
# chr9  = 255,204,0
# chr10 = 255,255,0
# chr11 = 204,255,0
# chr12 = 0,255,0
# chr13 = 53,128,0
# chr14 = 0,0,204
# chr15 = 102,153,255
# chr16 = 153,204,255
# chr17 = 0,255,255
# chr18 = 204,255,255
# chr19 = 153,0,204
# chr20 = 204,51,255
# chr21 = 204,153,255
# chr22 = 102,102,102
# chr23 = 153,153,153
# chrX  = 153,153,153
# chr24 = 204,204,204
# chrY  = 204,204,204
# chrM  = 204,204,153
# chr0  = 204,204,153
# chrUn = 121,204,61
# chrNA = 255,255,255
