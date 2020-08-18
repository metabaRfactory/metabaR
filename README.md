
<!-- README.md is generated from README.Rmd. Please edit that file -->
![MetabarF banner](metabaRF.png)

*MetabaR* is an R package enabling the import, handling and processing metabarcoding data that have been already processed with bioinformatic pipelines. It provides functions to identify and filter common molecular artifacts produced during the experimental workflow.

This package can be easily used in combination with others packages commonly used in ecology (*vegan*, *ade4*, *ape*, *picante*, etc.), and provides flexible graphic systems based on *ggplot2* to vizualize the data under both ecological and experimental perspectives.

More specifically, *MetabaR* provides:

-   Import functions of DNA metabarcoding data processed with common bioinformatic pipelines.    
-   Functions to manipulate the suite of data tables one usually deals with when working with DNA metabarcoding.    
-   Functions of data curation that are absent from the bioinformatics pipelines (e.g. detecting and tagging potential molecular artifacts such as contaminants, dysfynctional PCRs, etc.)
-   Functions for visualizing the data under both ecological (e.g. type of samples) and experimental (e.g. type of controls, distribution across the PCR plate design) perspectives.

*MetabaR* is developed on GitHub:

<https://github.com/metabaRfactory/metabaRffe>

Overall overview
----------------

![MetabarF over](metabaRF_overview.png)

Installation
------------

*MetabaR-F* can be installed from github with:

``` r
install.packages("devtools")
devtools::install_github("metabaRfactory/metabaRffe")
```

Dependencies:
- graphical R packages: *igraph*, *ggplot2* and *cowplot*
- data format-related R pacakges: *reshape2*, *seqinr*, *biomformat*
- ecological/statistical R packages: *vegan*, *ade4*

Example
-------

This is a basic example of use:

``` r
library(metabaR-F)
data(soil_euk)
summary_metabarlist(soil_euk)
```
