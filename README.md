# PEA: A Proteomics Expression Analysis, Application

## Table of Contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Imputation](#imputation)
* [DET Corrector](#detcorrector)
* [Protein Vector Normalization](#pvn)
* [Data Analytics](#dataanalytics)
<br><br>

## Introduction <a name="introduction"></a>
Written in Python and R - allows user to impute missing data using a variety of algorithms, and then estimate differential abundance on the corrected matrix. The proteomics matrix is computed in Python using the model of choice among missForest, KNN or RegImpute from the missingpy and StaticImpute libraries. All PCA and Volcano plots are generated in R using Shiny.
<br><br>

## Installation <a name="installation"></a>
Installation of Python3, R, and RStudio is required.
Python libraries can be installed via pip install -r requirements.txt.
R Libraries and App can be installed and loaded with the devtools github install command.
<br><br>

## Imputation <a name="imputation"></a>
By generating results from 3 different models, KNN K-nearest neighbors, missForest, and RegImpute, we are able to decide which model represents the data best. KNN, applied first, is an algorithm based on K-nearest neighbors, or the proximity between data points. In our case k=10, thus, there are 10 points forming clusters that are used to determine this proximity. Data points are computed based on the cluster that they belong to. MissForest is an algorithm based on iterations of predicting the missing values. Each time a prediction is made for a particular missing value, the algorithm uses this as a new training set to predict the next value for the same missing value. Data points are calculated based on previous predictions. RegImpute uses values to calculate missing values based on a linear regression of the current dataset. Missing values in the dataset are computed based on points within the linear regression that closely resemble the points in the neighboring technical replicates. While the new dataset cannot match the original dataset completely we are looking for a new matrix that introduces minimal change with high accuracy in terms of nearest technical replicates. While it is important to maintain the integrity of the original dataset, we plot the difference and note any outliers that may have changed. The results of some models may create false fold changes or p-values that are artificially inflated, the model creating these results are rejected. The model that produces changes in any PSM fold change aligned with the assigned protein fold change remains as a solid parameter for the PEA tool. Each dataset may be different, so the user can decide which model best represents their data.
<br><br>

## DET Corrector <a name="detcorrector"></a>
DET Corrector, corrects the peptide spectral matches based on Detection, Experimental and Technical error. Each abundance is normalized by the average of the fitting to of each model:
Technical Error, &lambda; = &gamma; x^3 + &alpha; x^2 + &beta; x + n
Experimental Error, &psi; = &gamma; x^3 + &alpha; x^2 + &beta; x + n
Detection Error, &omega; = &gamma; x^4 + &alpha; x^3 + &beta; x^2 + x + n
<br><br>

## Protein Vector Normalization <a name="pvn"></a>
PVN Normalizes the Protein Matrix, based on the Control and Treatment Variance; as well as, the Delta m/z between Proteins.
<br><br>

## Data Analytics <a name="dataanalytics"></a>
PEA starts out by generating a new matrix based on a percentage matrix and a chosen model. Then filters out all high confidence PSM matches with unique peptides at a number greater than two per protein. By uploading the new matrix into the Shiny app and choosing the channels for treatment and control, we can then move forward in the analyses. We analyze the noise with each channel through box-whisker plots. Each condition is then compared with a PCA plot, making sure that each channel within treatment and control does not contain any outliers. Volcano plots are produced based on significant fold changes and p-values. All plots and tables are saved for further evaluation.

After methods are applied for calculating and adjusting missing data in TMT Proteomics data, the file is further filtered with protein FDR confidence is high, unique peptides greater than 2, master proteins only, and no contaminants. Some of the graphs and tables produced include PCA plots, Volcano plots, and tables including all the statistics presented in the graphs. Applied here is a VSN normalization computed on the imputed matrix using a robust variant of the maximum-likelihood estimator for an additive-multiplicative error model and affine calibration. The model incorporates dependence of the variance on the mean intensity and a variance stabilizing data transformation. A linear model is fitted to the expression data for control and treatment, then t-statistics are computed by empirical Bayes moderation of standard errors towards a common value.
<br><br>