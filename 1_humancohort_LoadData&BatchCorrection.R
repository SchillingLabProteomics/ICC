# Load Data

library(tidyr)
library(readxl)
library(devtools)
library(diann)
library(tidyverse)

# Read and prepare DIA-NN output
dataraw <- diann_load("Data/ICC_remeasurement_TimsTOF_predLib_PSM.tsv")

#remove filepath from File.Name
dataraw <- mutate(dataraw, File.Name = str_replace(dataraw[[1]], ".*\\\\", ""))
dataraw <- mutate(dataraw, File.Name = str_replace(dataraw[[1]], ".raw.mzml$", ""))

#extract MQ-LFQ data from file - Genes.MQLFQ is already normalized via DIA-NN!

precursors <- diann_matrix(dataraw, pg.q = 0.01, proteotypic.only = T)
unique.genes <- diann_matrix(dataraw, id.header="Protein.Ids", 
                             quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01)
write.csv(unique.genes, "results_1_humancohort_LoadData//TimsTOF_ICC_humancohort_matrix.tsv")
