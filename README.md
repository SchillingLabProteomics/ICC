# Intrahepatic Cholangiocarcinoma (ICC) Proteomics Analysis Pipeline

## Overview

This repository contains the analysis pipeline for the manuscript:

"Proteomic Characterization of Intrahepatic Cholangiocarcinoma Identifies Risk-Stratifying Subgroups, Proteins Associated with Time-To-Recurrence, and mTOR Effector Molecule EIF4A1 as a Druggable Therapeutic Target."

The study is based on proteomic data from a human cohort of intrahepatic cholangiocarcinoma (ICC) comprising 80 tumor and 77 tumor-adjacent, non-malignant (TANM) tissue samples from 80 patients. The cohort was analyzed using DIA-NN 1.7 with a fully tryptic database, and a semi-specific database to study proteolytic events. Additionally, we analyzed 9 patient-derived ICC xenografts (PDX).


## How to Use the Code

The scripts are numbered to correspond to the main steps in the analysis pipeline described in the manuscript.

### 1. Data Preparation & Batch Correction

  `1_humancohort_LoadData&BatchCorrection.R`

  Loads DIA-NN output and performs batch correction.
  Outputs expression matrices for further analysis.

### 2. Unsupervised Analysis of Human Cohort

`2_humancohort_unsupervisedStats.Rmd`

Conducts PCA and hierarchical clustering.
Generates exploratory plots.

### 3. Tumor vs. TANM Comparison

`3_humancohort_Tumor-TANM.R`

Performs differential expression analysis between tumor and tumor-adjacent non-malignant tissues (TANM) using limma.

### 4. Unsupervised Analysis of Tumor Samples

`4_humancohort_tumor_unsupervisedStats.Rmd`

Clustering and PCA of tumor samples.

### 5. ECM vs. Turnover Cluster Analysis

`5_humancohort_ECM-TurnoverCluster.R`

Differential expression analysis between ECM-cluster and protein-turnover-cluster using limma.
Volcano plots and Kaplan-Meier survival analysis.

### 6. Survival Analysis

`6_humancohort_CPHM.Rmd`

Cox proportional hazards models (CPHM) to identify proteins associated with time-to-recurrence.

### 7. Semi-Specific Analysis

`7_humancohort_semi-specificAnalysis.Rmd`

Analysis of semi-specific peptides including abundance and PLS-DA.

### 8. Xenograft Analysis

`8_xenografts.Rmd`

Two-species analysis of proteomic data from patient-derived ICC xenografts (PDX).
PCA-based co-regulation analysis.


## Input Data

DIA-NN output files: Download from MassIVE Repository and place them in the Data/ folder. Reviewer login credentials for MassIVE are outlined in the Data Sharing section of the manuscript. The DIA-NN output files are: 



Expression Matrices: Available in Data/, for users who want to start with statistical analyses directly.

## Notes

The classifier model is not included in this repository as it is being prepared for a separate publication. It will be made publicly available before the final acceptance of this manuscript.
