# Intrahepatic Cholangiocarcinoma (ICC) Proteomics Analysis Pipeline

## Overview

This repository contains the analysis pipeline for the manuscript:

"Proteomic Characterization of Intrahepatic Cholangiocarcinoma Identifies Risk-Stratifying Subgroups, Proteins Associated with Time-To-Recurrence, and mTOR Effector Molecule EIF4A1 as a Druggable Therapeutic Target."

The study is based on proteomic data from a human cohort of intrahepatic cholangiocarcinoma (ICC) comprising 80 tumor and 77 tumor-adjacent, non-malignant (TANM) tissue samples from 80 patients. All samples were measured on a TimsTOF Flex coupled to an Evosep One chromatographic system. Data from earlier measurements on a Thermo QExactive coupled to an Easynlc nanoflow system is available but was not included in the publication. The cohort was analyzed using DIA-NN 1.9.2 with a fully tryptic database, and a semi-specific database to study proteolytic events. Additionally, we analyzed 9 patient-derived ICC xenografts (PDX).


## How to Use the Code

The scripts are numbered to correspond to the main steps in the analysis pipeline described in the manuscript.

### 1. Data Preparation & Batch Correction

  `1_humancohort_LoadData.R`

  Loads DIA-NN output and generate expression matrices for further analysis.

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

### 4B. Unsupervised Analysis of TANM Samples

`4B_humancohort_TANM_unsupervisedStats.Rmd`

Clustering and PCA of TANM samples.

### 5. ECM vs. Proliferation Cluster Analysis

`5_humancohort_ECM-ProliferationCluster.R`

Differential expression analysis between ECM-cluster and proliferation-cluster using limma.
Volcano plots and Kaplan-Meier survival analysis.

### 5B.TANM Cluster Analysis

`5B_humancohort_TANMCluster.R`

Differential expression analysis between TANM-clusters incl. Gene Set Enrichment.

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

All mass spectrometric raw data and analysis results are deposited in the MassIVE repository together with the relevant clinical annotation. Reviewer login credentials for MassIVE are outlined in the Data Sharing section of this manuscript. Once you are logged in, click "Browse Dataset Files". Please download the files outlined below and copy them into the data folder before executing the code. 

DIA-NN output files *(Filter by Collection = quant)*:<br>
`ICC_humancohort_DIA-NNoutput.tsv` `ICC_humancohort_semi-specific_DIA-NNoutput.tsv` <br> `ICC_xenografts_DIA-NNoutput.tsv`

Expression Matrices *(Filter by Collection = search)*:<br>
`TimsTOF_ICC_humancohort_matrix.tsv` `ICC_remeasurement_TimsTOF_semispec_reannot.tsv` <br> `semi-specific_precursors_sequence+ID+Abundance.csv` `ICC_xenografts_matrix.csv` `ICC_xenografts_matrix_humanproteins.tsv`

Analysis-specific annotation *(Filter by Collection = metadata)*:<br>
`annotation_batch_correction.xlsx` `study_design.xlsx`

Charts and Tables *(Filter by Collection = metadata; also provided as supplementary tables during manuscript submission)*:  
`Table S1 - Patient Table.xlsx` `Table S2 - Coxph_adj.hits_stats.xlsx` `Table S3 - Cox_stats_table.xlsx` <br> `Table S4 - PDX_treatment_eft226.xlsx` `Table S5 - ICC_proteomic_data_after_batch_correction.xlsx` <br> `Table S6 - PDX proteomic data.xlsx` `Table S7 - Differential Proteins Clusters.xlsx` <br> `Table S8 - Differential Proteins Tumor TANM .xlsx` 

## Notes

The classifier model is not included in this repository as it is being prepared for a separate publication.
