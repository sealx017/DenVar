DenVar
================

This is a R package for visualization and supervised cell annotation of
Vectra Polaris and MIBI data. The package provides basic marker
visualization tools like: heat-map and ridge plot. Given a set of
annotated training images, it can also perform a random forest based
cell phenotyping for non-annotated images.

## Loading required packages

First, we load our package: DenVar and a few other required packages.
One can install the developmental version of VectraMIBI by running the
command: **devtools::install\_github(‘sealx017/DenVar’)**.

``` r

require(DenVar)
# Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
# logical.return = TRUE, : there is no package called 'DenVar'
require(pheatmap)
require(ggplot2)
require(LaplacesDemon)
```

## Loading the dataset

Next, we import the training dataset named as “train\_data”, which has
marker data (34 many) for different cells from two different images. The
first two columns: “SampleID” and “cellLabelInImage” of the training
dataset respectively correspond to the image number and the cell label
inside that image. The third column (“Group”) corresponds to the
annotated cell type.

``` r
PD1_data = read.csv("Data/PD1_TNBC_MIBI.csv")
clinical_data = read.csv("Data/Survival_TNBC_MIBI.csv")
```

## Computing the KDEs and the JSD matrix

Next, we import the training dataset named as “train\_data”, which has
marker data (34 many) for different cells from two different images. The
first two columns: “SampleID” and “cellLabelInImage” of the training
dataset respectively correspond to the image number and the cell label
inside that image. The third column (“Group”) corresponds to the
annotated cell type.
