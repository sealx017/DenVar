DenVar ![Image](Denvar_logo.tif)
================

This is a R package implementing the proposed method from the paper
“DenVar: Density-based Variation analysis of multiplex imaging data”.

## Loading required packages

First, we load our package, DenVar and a few other required packages.
One can install the developmental version of DenVar by running the
command: **devtools::install\_github(‘sealx017/DenVar’)**.

``` r
#devtools::install_github('sealx017/DenVar')
require(DenVar)
require(LaplacesDemon)
require(pheatmap)
require(ggplot2)
require(survival)
require(survminer)
require(coxme)
```

## Loading the dataset

Next, we import the example datasets named, “PD1\_TNBC\_MIBI.csv” and
“Survival\_TNBC\_MIBI.csv” which are respectively the PD-1 marker
expression data and the survival data of 33 patients from TMBC MIBI
data. Both the datasets have a column named “SampleID” denoting
individual patient IDs.

``` r

PD1_data = read.csv("Data/PD1_TNBC_MIBI.csv")[,-1]
clinical_data = read.csv("Data/Survival_TNBC_MIBI.csv")[,-1]
```

## Computing the KDEs and the JSD matrix

We compute the KDE of PD-1 expression in all the patients and store it
in an array. We keep number of grid points (ngrids) in the KDE at 1024.
Using the array of KDEs, we compute pairwise Jensen-Shannon Distance
(JSD) between the patient densities. We create the JSD matrix having
distances between all possible patient pairs.

``` r

KDE_of_all = Array_KDE(Data = PD1_data, ngrids = 1024)
JSD_matrix_between_all = JSD_matrix(Array_dens = KDE_of_all, Data = PD1_data)
```

## Heatmap of the JSD matrix

Next, we visualize the computed JSD matrix as a heatmap.

``` r

print(pheatmap(JSD_matrix_between_all, cluster_rows = TRUE,
cluster_cols = TRUE,show_colnames = T, show_rownames = T,treeheight_row=0,treeheight_col = 0))
```

<img src="README_files/figure-gfm/heatmap of JSD matrix-1.png" width="100%" />

## Hierarchical clustering based on the JSD matrix

We use the JSD matrix for hierarchical clusetring of the patients into
two meanigful clusters.

``` r

hier_clustering = hclust(as.dist(JSD_matrix_between_all), method="ward.D")
hier_groups = as.data.frame(cutree(hier_clustering, k=2)); colnames(hier_groups) = "Cluster"
```

## CoxPH model and KM plots

Using the vector of cluster-labels, we fit a CoxPH model and plot the
Kaplan-Meyer (KM) curves stratified by the clusters. We display the
estimated Hazard Ratio (HR) and the corresponding p-value as well.

``` r

surv = cbind(clinical_data, hier_groups)
coxplot = coxPH_plot(surv)
```

<img src="README_files/figure-gfm/CoxPh model with the clusters-1.png" width="100%" />

``` r
#print(coxplot$plot)
```

## CoxPH model with random effects (Coxme)

We use the computed JSD matrix directly in a CoxPH model with random
effect (also known as Coxme). We show the estimated variance component
and the corresponding LRT for testing its significance.

``` r

res.coxme = coxme_model(surv, JSD_matrix_between_all)
print(res.coxme)
# $Variance_estimate
# $Variance_estimate$SampleID
#   Vmat.1 
# 2.239131 
# 
# 
# $LRT
# Integrated 
#   3.603008 
# 
# $LRT_pvalue
# Integrated 
# 0.05767512
```
