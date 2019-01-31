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

## A few example plots created with karyoploteR

These images are all created with karyoploteR and are part of the documented 
examples in the [karyoploteR's tutorial and examples page](https://bernatgel.github.io/karyoploter_tutorial/).
Click on them to see how the code needed to create them.


<p align="center">
  <a href="https://bernatgel.github.io/karyoploter_tutorial/Examples/GeneExpression/GeneExpression.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial/Examples/GeneExpression/images/Figure13-1.png" width="40%"  alt="A karyoploteR example plotting differential expression results computed with RNA-seq data from Drosophila Melanogaster" title="Differential expression results computed with RNA-seq data from Drosophila Melanogaster" style="max-width:100%;margin-right:5%;" ></img>
  </a>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneDensityIdeograms/GeneDensityIdeograms.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneDensityIdeograms/images/Figure9-1.png" width="40%"  alt="A karyoploteR example plotting the density of genes instead of the ideograms" title="Usiong the density of genes instead of ideograms"></img> 
  </a>
  <br>
  <br>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/NucleotideFrequency/NucleotideFrequency.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/NucleotideFrequency/images/Figure11-1.png" width="40%"  alt="A karyoploteR example plotting the nucleotide frequency, genes and CpG-islands on a small genomic region" title="The nucleotide frequency, genes and CpG-islands on a small genomic region"></img>
  </a>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/SNPArray/SNPArray.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/SNPArray/images/Figure4-1.png" width="40%"  alt="A karyoploteR example plotting raw SNP-array data" title="Raw SNP-array data"></img>
  </a>
  <br>
  <br>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/PVivaxGenes/PVivaxGenes.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/PVivaxGenes/images/Figure6-1.png" width="40%"  alt="A karyoploteR example plotting the genes from Plasmodium Vivax PvP01 genome version" title="The genes from Plasmodium Vivax PvP01 genome version"></img>
  </a>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/PlotGenes/PlotGenes.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/PlotGenes/images/Figure2-1.png" width="40%"  alt="A karyoploteR example plotting genes positioned on the genome" title="Genes positioned on the genome"></img>
  </a>
  <br>
  <br>
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/Rainfall/Rainfall.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/Rainfall/images/Figure3-1.png" width="40%"  alt="A karyoploteR example plotting a rainfall plot showing the distances between consecutive somatic variants" title="A rainfall plot showing the distances between consecutive somatic variants"></img>
  </a>  
  <a href="https://bernatgel.github.io/karyoploter_tutorial//Examples/CpGIslands/CpGIslands.html" target="_blank">
    <img src="https://bernatgel.github.io/karyoploter_tutorial//Examples/CpGIslands/images/Figure5-1.png" width="40%"  alt="A karyoploteR example plotting the density and positions of CpG islands along the genome" title="The density and positions of CpG islands along the genome"></img>
  </a>
</p>












