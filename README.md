DenVar
================

This is a R package implementing the proposed method from the paper
“DenVar: Density-based Variation analysis of multiplex imaging data”.

## Loading required packages

First, we load our package, DenVar and a few other required packages.
One can install the developmental version of DenVar by running the
command: **devtools::install_github(‘sealx017/DenVar’)**.

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

We import the example datasets named, “Marker_data_TNBC_MIBI.csv” and
“clinical_data_TNBC_MIBI.csv” extracted from the TNBC data used in Keren
et al. (2018). The former has the PD1 marker-expression data in 170, 171
cells from 33 subjects. The latter has recurrence and survival data of
those same subjects. Both the datasets have a column named “ID” denoting
the individual subject IDs.

``` r

PD1_data = read.csv("Data/Marker_data_TNBC_MIBI.csv")
clinical_data = read.csv("Data/clinical_data_TNBC_MIBI.csv")
```

## Computing the KDEs and the JSD matrix

We compute the KDE of PD-1 expression in all the subjects and store it
in an array. We keep number of grid points (ngrids) in the KDE at 1024.
Using the array of KDEs, we compute the pairwise Jensen-Shannon Distance
(JSD) between the subject-specific densities. We create the JSD matrix
having distances between all possible patient pairs.

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

We subject the JSD matrix to hierarchical clustering to group the
subjects into two meaningful clusters.

``` r

hier_clustering = hclust(as.dist(JSD_matrix_between_all), method="ward.D")
hier_groups = as.data.frame(cutree(hier_clustering, k=2)); colnames(hier_groups) = "Cluster"
```

## CoxPH model with Recurrence Outcome and KM plot

Using the recurrence outcome (time and censoring-indicator), we fit a
CoxPH model and plot the Kaplan-Meyer (KM) curves stratified by the
clusters. We display the estimated Hazard Ratio (HR) and the
corresponding p-value as well. Note that we have created a matrix named
“surv” by merging the outcome information and the vector of
cluster-labels and changed the names of the censoring-indicator and time
columns to “Censored” and “Time” respectively.

``` r

surv = cbind(clinical_data[,c(1:3)], hier_groups)
colnames(surv)[2:3] = c("Censored", "Time")
print(head(surv))
#   ID Censored Time Cluster
# 1  1        1    9       1
# 2  2        0  745       1
# 3  3        0 3130       2
# 4  4        1   31       2
# 5  5        0 1683       2
# 6  6        0 2275       1
coxplot = coxPH_plot(surv, name = "Non-recurrence probability")
```

<img src="README_files/figure-gfm/CoxPh model with the clusters-1.png" width="100%" />

``` r
#print(coxplot$plot)
```

## CoxPH model with random effects (Coxme) with Recurrence Outcome

Using the recurrence outcome, we use the computed JSD matrix directly in
a CoxPH model with random effect (also known as Coxme). We show the
estimated variance component and the corresponding LRT for testing its
significance.

``` r

res.coxme = coxme_model(surv, JSD_matrix_between_all)
print(res.coxme)
# $Variance_estimate
# $Variance_estimate$ID
#   Vmat.1 
# 0.658883 
# 
# 
# $LRT
# Integrated 
#   1.130152 
# 
# $LRT_pvalue
# Integrated 
#   0.287743
```

## CoxPH model with Survival Outcome and KM plot

Using the survival outcome (time and censoring indicator), we fit a
CoxPH model and plot the Kaplan-Meyer (KM) curves stratified by the
clusters. We display the estimated Hazard Ratio (HR) and the
corresponding p-value as well. Note that we have created a matrix named
“surv” by merging the outcome information and the vector of
cluster-labels and changed the names of the censoring-indicator and time
columns to “Censored” and “Time” respectively.

``` r

surv = cbind(clinical_data[,c(1, 4:5)], hier_groups)
colnames(surv)[2:3] = c("Censored", "Time")
print(head(surv))
#   ID Censored Time Cluster
# 1  1        1 2612       1
# 2  2        1  745       1
# 3  3        0 3130       2
# 4  4        0 2523       2
# 5  5        0 1683       2
# 6  6        0 2275       1
coxplot = coxPH_plot(surv, name = "Survival probability")
```

<img src="README_files/figure-gfm/CoxPh model with the clusters and survival-1.png" width="100%" />

``` r
#print(coxplot$plot)
```

## CoxPH model with random effects (Coxme) with Survival Outcome

Using the survival outcome, we use the computed JSD matrix directly in a
CoxPH model with random effect (also known as Coxme). We show the
estimated variance component and the corresponding LRT for testing its
significance.

``` r

res.coxme = coxme_model(surv, JSD_matrix_between_all)
print(res.coxme)
# $Variance_estimate
# $Variance_estimate$ID
#    Vmat.1 
# 0.4727335 
# 
# 
# $LRT
# Integrated 
#   0.318442 
# 
# $LRT_pvalue
# Integrated 
#  0.5725455
```
