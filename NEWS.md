# karyoploteR 1.7.1

* Added `kpPlotBigWig` to plot data in bigwig files, usually data derived from bam coverage for ChIP-seq, etc...


## Other

* Added a `digits` parameter to `kpAddBAseNumbers` to control the number of digits after the decimal point in genome position lables


# karyoploteR 1.5.4

## New features

* Added `kpPlotGenes` and `kpPlotTranscripts` to plot gene and transcript models
* Added `kpArea` to plot shaded areas. Ideal for coverage plots, RNA-seq, ChIP-seq, etc...

## API changes


## Other

* Added plot.type=4 to `plotDefaultPlotParams`
* Added new examples and extended the tutorial at the karyoploteR tutorial and examples site at https://bernatgel.github.io/karyoploter_tutorial/


## Bug Fixes

* Various minor bug fixes.
* Various documentation fixes



# karyoploteR 1.3.11

## New features

* Added zooming to create plots of regions smaller than a whole chromosome. 
* Added `kpPlotLoess` to plot a fitted loess and confidence interval for data points.
* Added `kpPlotRainfall` to create rainfall plots from variants.
* Added `kpPlotLinks` to plot connections between genome regions even in different chromosomes.


## API changes

* New default in `plotKaryotype`: now `plot.type` defaults to 1, a ideogram with a single data panel above it


## Other

* Added unit testing
* Created a karyoploteR tutorial and examples site at https://bernatgel.github.io/karyoploter_tutorial/


## Bug Fixes

* Fixed a bug causing a misalignment between data points and plotting parameters in some edge cases.
* Various minor bug fixes.

