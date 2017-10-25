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

