#INTERNAL

#Returns a table with the colors associated with each cytoband.
# By default it returns the colors defined in Circos.

getCytobandColors <- function(color.table) {
  
  if(is.null(color.table)) {
    color.table <- list(gneg="#FFFFFF",
                     gpos25="#C8C8C8",
                     gpos33="#D2D2D2",
                     gpos50="#C8C8C8",
                     gpos66="#A0A0A0",
                     gpos75="#828282",
                     gpos100="#000000",
                     gpos="#000000",
                     stalk="#647FA4", #repetitive areas
                     acen="#D92F27", #centromeres
                     gvar="#DCDCDC")
  }
  return(color.table)
}