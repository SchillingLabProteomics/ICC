### Script for running multi-group limma analysis ---- 
## Miguel Cosenza 02.07.2020 

### 1. Please provide a meaningful name for the dataset/experiment ----

exper_code <- "Limma TANM vs. tumor"

### 2. How many top significant hits do you want to plot (boxplots group vs intensity)? 

n_top_hits <- 12

### 3. Which groups do you want to compare?  

## Required packages ----

packages <- c("dplyr", "here", "tidyr", "ggplot2", "rmarkdown", "knitr", "reshape")

biopackgs <- c("limma")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
          install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
          if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
          
          BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
          
}

library(dplyr)
library(stringr)
library(limma)
library(rmarkdown)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(reshape)
library(plotly)
library(ggrepel)
library(readxl)
library(tidyverse)
`%notin%` <- Negate(`%in%`)

## Load data ----
#modify annot_dat processing to change comparison

annot_dat <- read_excel("Data/Table S1 - Patient Table.xlsx") %>%
  rename(Sample_ID = 'Sample ID') %>%
  filter(Sample_ID %notin% c("KB106", "KB42", "KB147", "KB53", "KB213", "KB215", "KB55", "KB23", "KB209")) 

expr_dat <- read_excel("Data/Table S5 - ICC_proteomic_data_after_batch_correction.xlsx", na = "NA") %>%
            dplyr::rename(ID = `...1`) %>%
            dplyr::select(ID, one_of(annot_dat$Sample_ID))

# Filter for 20% missing values
expr_dat <- filter(expr_dat, rowSums(!is.na(expr_dat)) / (ncol(expr_dat)-1) > 0.8)

# Median Normalization
expr_dat <- mutate_if(expr_dat, is.numeric, function(x, na.rm = TRUE) x - median(x, na.rm = TRUE))

# Load human proteome
hupo <- read.csv("Data/uniprot-proteome UP000005640+reviewed yes(incl.GO).tab", sep = "\t") %>%
  rename(Protein = "Entry")
entries <- str_split_fixed(hupo$Entry.name, "_", 2)
hupo$Entry.name <- entries

## Define the design matrix ----

groups <- as.factor(annot_dat$Tissue)

design <- model.matrix(~0+groups)

row.names(design) <- annot_dat$Sample_ID

colnames(design) <- colnames(design) %>% str_remove(., "groups") %>% str_trim()

# 3. DEFINE WHICH GROUPS YOU WISH TO COMPARE ----

contrast.matrix <- makeContrasts(Tumor - normal, levels=design) 

## Prep expression data matrix ----
n_contrasts <- dim(contrast.matrix)[2]

tomat <- dplyr::select(expr_dat,
                       -ID, row.names(design)) %>% as.data.frame()

row.names(tomat) <- expr_dat$ID

mat <- as.matrix(tomat)

## Execute the linear model / limma ----

fit <- lmFit(mat, design = design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

output_limma2 <- topTable(fit2, adjust.method = "BH", number = Inf)

output_limma2$Protein <- row.names(output_limma2)


## Generate output ----

## output tabular list ----

list_tabular <- list()

for (i in 1:n_contrasts){
          
          outlim <- topTable(fit2, coef = i, adjust.method = "BH", number = Inf)
          outlim$Protein <- row.names(outlim)
          outlim$Contrast <- colnames(contrast.matrix)[i]
          list_tabular[[i]] <- outlim
}

names(list_tabular) <- colnames(contrast.matrix)

if(!dir.exists("Output")){dir.create("Output")}

write.table(output_limma2,
            file = "Output/tab_output_general_Ftest_limma.txt",
            row.names = FALSE, col.names = TRUE)

sig_hits <- dplyr::filter(output_limma2, #MODIFIED!!! 
                               adj.P.Val <= 0.01) %>% row.names(.)

n_significant <- length(sig_hits)

top_n_hits <- slice_min(output_limma2, n = n_top_hits, order_by = adj.P.Val, with_ties = FALSE)


## Extract basic information about each contrasts ----
getinfo <- function(x){
          x %>% dplyr::filter(adj.P.Val <= 0.01) %>% dim(.) %>% .[1]
}

## Prep some boxplots for the top proteins with lowest P-values ----

slim_expr <- pivot_longer(expr_dat, cols = colnames(mat),
                          names_to = "Sample_ID",
                          values_to = "Abundance") 

slim_expr_g <- left_join(slim_expr, annot_dat,
                         by = "Sample_ID")

hits_expr <- filter(slim_expr_g,
                    ID %in% row.names(top_n_hits))


boxplots <- ggplot(hits_expr,
                   aes(x = Tissue, y = Abundance)) +
          geom_boxplot()+
          facet_wrap(ID ~ .)

print(boxplots)

## Prep volcano plots for contrasts ----

merged_limmacontr <- reshape::merge_all(list_tabular)

write.table(merged_limmacontr,
            file = "Output/tab_output_per_contrasts_limma.txt",
            row.names = FALSE, col.names = TRUE)


tovolc <- merged_limmacontr %>% 
          mutate(Differentially_expressed = case_when(adj.P.Val <= 0.01 ~ TRUE, TRUE ~ FALSE)) #MODIFIED
tovolc <- left_join(tovolc, hupo, by = "Protein") 

tophits <- slice_max(tovolc, order_by = logFC, n = 10) %>%
  bind_rows(slice_min(tovolc, order_by = logFC, n = 10))

write.csv(tovolc, file = "Output/Limmaresult_inclUNIPROTannotation.csv")

#simple volcano
# simple volcano plot of patient matched analysis 
tovolc$Gene.names <- str_split_fixed(tovolc$Gene.names, " ", n = Inf)
tovolc$Gene.names <- tovolc$Gene.names[,1]

simpvol2 <- ggplot(data = tovolc, mapping = aes(x = logFC, y = -log10(adj.P.Val), label = Gene.names)) +
  geom_point(data = filter(tovolc,  adj.P.Val > 0.05),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "black", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC > 0, adj.P.Val <= 0.05),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "red", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC < 0, adj.P.Val <= 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "blue", alpha = 0.4)+
  geom_text_repel(aes(label = ifelse(Protein %in% tophits$Protein, as.character(Gene.names),'')),
                  hjust = 0, vjust = 0, size = 3, max.overlaps = Inf,
                  min.segment.length = 3, point.size = NA, show.legend = FALSE) +
  theme_minimal() +
  #labs(title = "Tumor / TANM") +
  ylab("-log10 (adj. p-value)") +
  xlab("log2 fold change")

print(simpvol2)
