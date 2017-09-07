[![Build Status](https://travis-ci.org/bernatgel/karyoploteR.svg?branch=master)](https://travis-ci.org/bernatgel/karyoploteR)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

# karyoploteR - An R/Biocondutor package to plot arbitrary data along the genome


![karyoploteR Example](https://raw.githubusercontent.com/bernatgel/karyoploter_examples/master/Examples/Examples/MultipleDataTypes/figure/Figure-1.png "Example of plot created with karyoploteR")

## Description

[karyoploteR](http://bioconductor.org/packages/karyoploteR) is an R package to 
plot data along the genome using a karyotype style plot.

It is entirely based on R base graphics and inspired by the R base graphics API. 
It includes functions to plot primitive graphic elements such as points, lines,
rectangles, text, etc mapped into the genome plot coordinatesand and higher 
level functions to plot a heatmap, the regions in a GenomicRanges object
or the cumulative coverage of such regions.

Data positioning and track configuration has been inspired by Circos and does
not explicitly understands the concept of track. Thus, it is possible to freely specify 
where to plot the data and to create plots with multiple independent tracks or
overlapping representations.

It is highly configurable and in addition to the parametrizatiopn of the 
different data plotting functions, it is possible to specify custom functions 
for every plotting action from the basic chromosome bands to the chromosome labels
or base numbers as well as creating completely new plotting functions.

## How to use it

Documentation ([vignette](http://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.pdf) and [user manual](http://bioconductor.org/packages/devel/bioc/manuals/karyoploteR/man/karyoploteR.pdf)) is available at the karyoploteR's 
Bioconductor landing page at [http://bioconductor.org/packages/karyoploteR](http://bioconductor.org/packages/karyoploteR)

## Tutorial and Examples

In addition to the documentation above, a short tutorial and some examples can be found at [https://bernatgel.github.io/karyoploter_tutorial/](https://bernatgel.github.io/karyoploter_tutorial/)

## <a name="Citing"></a>Citing karyoploteR

karyoploteR has been developed by [Bernat Gel](https://twitter.com/bernatgel) and [Eduard Serra](mailto:eserra@igtp.cat) at [IGTP](http://www.germanstrias.org/)
Hereditary Cancer Group.

If you use karyoploteR in your research, please cite the [Bioinformatics paper](https://academic.oup.com/bioinformatics/article/3857734/karyoploteR-an-R-Bioconductor-package-to-plot) describing it:

Bernat Gel & Eduard Serra. (2017). *karyoploteR: an R/Bioconductor package to plot customizable genomes displaying arbitrary data*. Bioinformatics, 31–33. [doi:10.1093/bioinformatics/btx346](https://doi.org/10.1093/bioinformatics/btx346)





