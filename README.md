# karyoploteR - An R/Biocondutor package to plot arbitrary data along the genome

## Description

[karyoploteR](http://bioconductor.org/packages/karyoploteR) is an R package to 
plot data along the genome using a karyotype style plot.

It is entirely based on R base graphics and inspired by the R base graphics API. 
It includes functions to plot primitive graphic elements such as points, lines,
rectangles, text, etc mapped into the genome plot coordinatesand and higher 
level functions to plot a heatmap, the regions in a GenomicRanges object
or the cumulative coverage of such regions.

Data positioning and track configuration has been inspired by Circos and does
not explicitely understands the concept of track. Thus, it is possible to freely specify 
where to plot the data and to create plots with multiple independent tracks or
overlapping representations.

It is highly configurable and in addition to the parametrizatiopn of the 
different data plotting functions, it is possible to specify custom functions 
for every plotting action from the basic chromosome bands to the chromosome labels
or base numbers as well as creating completely new plotting functions.

## How to use it

Documentation ([vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf) and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf)) is available at the karyoploteR's 
Bioconductor landing page at [http://bioconductor.org/packages/karyoploteR](http://bioconductor.org/packages/karyoploteR)

## Examples

Some usage examples can be found at [https://github.com/bernatgel/karyoploter_examples](https://github.com/bernatgel/karyoploter_examples)

![karyoploteR Example](https://raw.githubusercontent.com/bernatgel/karyoploter_examples/master/Examples/CompleteExamples/MultipleDataTypes/figure/Figure-1.png "Example of plot created with karyoploteR")

