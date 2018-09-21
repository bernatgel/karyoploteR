# karyoploteR 1.7.9

## Bug Fixes

* Fixed a bug in the coordinate change function where plotting was out of place
or even invisible if the zoom object had addition seqlevels

# karyoploteR 1.7.4

## Other

* New parameter in kpPlotMarkers to allow labels to move beyond the 
chromosome limits when repositioning to avoid label overlaps



# karyoploteR 1.7.4

## Other

* The zoom region in plotKaryotype can be specified in any format accepted
by regioneR::toGRanges, including UCSC/IGV style "chr9:23000-40000".


# karyoploteR 1.7.3

* Added `kpPlotBAMCoverage` to plot the exact coverage from a BAM file

## Other

* Improved performance of kpPlotBAMDensity. Specially in zoomed plots.

## Bug Fixes

* kpAxis: Axis were not visible in zoomed plots. They are now visible.


# karyoploteR 1.7.1

* Added `kpPlotBigWig` to plot data in bigwig files, usually data derived from BAM coverage for ChIP-seq, etc...


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

