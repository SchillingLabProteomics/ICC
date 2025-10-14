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
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
`%notin%` <- Negate(`%in%`)

## Load data ----

annot_dat <- read_excel("Data/Table S1 - Patient Table.xlsx") %>%
  dplyr::rename(Sample_ID = 'Sample') %>%
  filter(Sample_ID %notin% c("KB106", "KB42", "KB147", "KB53", "KB213", "KB215", "KB55", "KB23", "KB209")) 

expr_dat <- read.csv("results_2_humancohort_unsupervisedStats/TimsTOF_ICC_ExprMat_log2_median.csv") %>%
  dplyr::select(-X) %>%
  dplyr::select("Protein", one_of(annot_dat$Sample_ID))

annot_dat <- filter(annot_dat, Sample_ID %in% colnames(expr_dat))

# Filter for 20% missing values
expr_dat <- filter(expr_dat, rowSums(!is.na(expr_dat)) / (ncol(expr_dat)-1) > 0.8) 

# Load human proteome
hupo <- read.csv("Data/uniprot-proteome UP000005640+reviewed yes(incl.GO).tab", sep = "\t") %>%
  dplyr::rename(Protein = "Entry")
entries <- str_split_fixed(hupo$Entry.name, "_", 2)
hupo$Entry.name <- entries

## Define the design matrix ----

groups <- as.factor(annot_dat$Group)

design <- model.matrix(~0+groups)

row.names(design) <- annot_dat$Sample_ID

colnames(design) <- colnames(design) %>% str_remove(., "groups") %>% str_trim()

# 3. DEFINE WHICH GROUPS YOU WISH TO COMPARE ----

contrast.matrix <- makeContrasts(Tumor - normal, levels=design) 

## Prep expression data matrix ----
n_contrasts <- dim(contrast.matrix)[2]

tomat <- dplyr::select(expr_dat,
                       -Protein, row.names(design)) %>% as.data.frame()

row.names(tomat) <- expr_dat$Protein

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

write.table(output_limma2,
            file = "results_3_humancohort_Tumor-TANM/tab_output_general_Ftest_limma.txt",
            row.names = FALSE, col.names = TRUE)

sig_hits <- dplyr::filter(output_limma2, 
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
                    Protein %in% row.names(top_n_hits))


boxplots <- ggplot(hits_expr,
                   aes(x = Group, y = Abundance)) +
          geom_boxplot()+
          facet_wrap(Protein ~ .)

print(boxplots)

## Prep volcano plots for contrasts ----

merged_limmacontr <- reshape::merge_all(list_tabular)

write.table(merged_limmacontr,
            file = "results_3_humancohort_Tumor-TANM/tab_output_per_contrasts_limma.txt",
            row.names = FALSE, col.names = TRUE)


tovolc <- merged_limmacontr %>% 
          mutate(Differentially_expressed = case_when(adj.P.Val <= 0.01 ~ TRUE, TRUE ~ FALSE)) #MODIFIED
tovolc <- left_join(tovolc, hupo, by = "Protein") 

tophits <- slice_max(tovolc, order_by = logFC, n = 10) %>%
  bind_rows(slice_min(tovolc, order_by = logFC, n = 10))

write.csv(tovolc, file = "results_3_humancohort_Tumor-TANM/Limmaresult_inclUNIPROTannotation.csv")

#simple volcano
# simple volcano plot of patient matched analysis 
tovolc$Gene.names <- str_split_fixed(tovolc$Gene.names, " ", n = Inf)
tovolc$Gene.names <- tovolc$Gene.names[,1]

simpvol2 <- ggplot(data = tovolc, mapping = aes(x = logFC, y = -log10(adj.P.Val), label = Gene.names)) +
  geom_point(data = filter(tovolc,  adj.P.Val > 0.01),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "black", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC > 0, adj.P.Val <= 0.01),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "#A52A2A", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC < 0, adj.P.Val <= 0.01), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "#4682B4", alpha = 0.4)+
  geom_text_repel(aes(label = ifelse(Protein %in% tophits$Protein, as.character(Gene.names),'')),
                  hjust = 0, vjust = 0, size = 3, max.overlaps = Inf,
                  min.segment.length = 3, point.size = NA, show.legend = FALSE) +
  theme_classic() +
  #labs(title = "Tumor / TANM") +
  ylab("-log10 (adj. p-value)") +
  xlab("log2 fold change")

print(simpvol2)

ggsave("results_3_humancohort_Tumor-TANM/volcano.pdf", plot = simpvol2, width = 63, height = 42, 
       units = "mm", dpi = 300, scale = 2)

tiff("results_3_humancohort_Tumor-TANM/volcano.tiff", units = "in", width = 10, height = 6, res = 300, pointsize =(10*100/72))
print(simpvol2)
dev.off()

geneList <- tovolc$logFC
names(geneList) = as.character(tovolc$Protein)
geneList = sort(geneList, decreasing = TRUE)

# gse <- readRDS("results_3_humancohort_Tumor-TANM/GeneSetEnrichment.RData")
gse <- gseGO(geneList, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "UNIPROT", minGSSize = 80, maxGSSize = 500,
             pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)

saveRDS(gse, file = "results_3_humancohort_Tumor-TANM/GeneSetEnrichment.RData")

ridgeplot <- ridgeplot(gse, label_format = 20) + 
  theme(axis.text.y = element_text(size = 6)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  labs(x = "log2-fold change", title = "Ontology: Biological Processes")

tiff("results_3_humancohort_Tumor-TANM/ridgeplot.tiff", units = "in", width = 8, height = 6, res = 300, pointsize =(8*100/72))
ridgeplot
dev.off()

plot <- gse@result 
plot1data <- slice_max(gse@result, enrichmentScore, n = 10) %>%
  bind_rows(slice_min(gse@result, enrichmentScore, n = 10)) %>%
  mutate(Description = str_wrap(Description, width = 50))

plot1 <- ggplot(plot1data, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore), size = setSize, color = p.adjust)) +
  geom_segment(aes(x = 0, xend = enrichmentScore, 
                   y = reorder(Description, enrichmentScore), yend = reorder(Description, enrichmentScore)),
               color = "grey", linewidth = 0.2) +
  geom_point() +
  theme_classic() +
  theme(axis.text.y= element_text(size = 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  geom_vline(xintercept = 0) +
  ylab("")
plot1

tiff("results_3_humancohort_Tumor-TANM/dotplot.tiff", units = "in", width = 10, height = 4, res = 300, pointsize =(8*100/72))
plot1
dev.off()

Tumor <- filter(tovolc, logFC >= 0)
geneListTumor <- Tumor$Protein
geneListTumor <- Tumor$logFC
names(geneListTumor) = as.character(Tumor$Protein)
geneListTumor = sort(geneListTumor, decreasing = TRUE)

TANM <- filter(tovolc, logFC <= 0)
geneListTANM <- TANM$Protein
geneListTANM <- TANM$logFC * -1
names(geneListTANM) = as.character(TANM$Protein)
geneListTANM = sort(geneListTANM, decreasing = TRUE)

bla <- list(TANM = geneListTANM, Tumor = geneListTumor)

data(gcSample)
ck <- compareCluster(bla, fun = gseGO, OrgDb = org.Hs.eg.db, keyType = "UNIPROT", minGSSize = 80)
plot2 <- dotplot(ck)
plot2

tiff("results_3_humancohort_Tumor-TANM/dotplot2.tiff", units = "in", width = 12, height = 4, res = 300, pointsize =(8*100/72))
plot2
dev.off()


# KEGG pathway enrichment
tovolcsign <- filter(tovolc, adj.P.Val <= 0.05)
input_pos <-  bitr(tovolcsign$Protein, fromType = "UNIPROT", 
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db) %>%
  dplyr::rename(Protein = "UNIPROT")

fgdata <- left_join(input_pos, dplyr::select(tovolcsign, "Protein", "logFC"), by = "Protein") %>%
  add_count(ENTREZID) %>%
  filter(n == 1) %>%
  dplyr::select(-n) %>%
  column_to_rownames(var = "ENTREZID") %>%
  dplyr::select(-Protein)

KEGG_enrichment_result <- enrichKEGG(input_pos$ENTREZID,
                                     organism = 'hsa',
                                     pvalueCutoff = 0.05)

hsa <- KEGG_enrichment_result$ID

KEGG_MN <- cbind(KEGG_enrichment_result$ID, KEGG_enrichment_result$Description, 
                 KEGG_enrichment_result$p.adjust, KEGG_enrichment_result$Count, 
                 KEGG_enrichment_result$geneID)

names(KEGG_MN) <- c("KeggPathway","Description","p-value","Count","Genes")
write.table(KEGG_MN,
            file = "results_3_humancohort_Tumor-TANM/KEGG_90_NM.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t"
)

setwd("results_3_humancohort_Tumor-TANM/KEGG_pathway/")
for(varhsa in hsa){
  pv.out <- pathview(gene.data =fgdata, 
                     pathway.id = varhsa, 
                     species = "hsa", 
                     out.suffix = "90perc_NM", 
                     kegg.native = TRUE,
                     same.layer = FALSE,
                     low = list(gene = "blue"), 
                     mid =list(gene = "gray"), 
                     high = list(gene = "red"))
}
setwd("D:/data/Klara/cohort_TimsTOF/analysis")

# Proteases
hupo2 <- read.delim("Data/uniprot-proteome UP000005640+reviewed yes.tab") %>%
  separate_rows(Gene.names, sep = " ")
proteases <- read_excel("Data/MEROPS_human_proteases.xlsx") %>%
  filter(!is.na(Gene)) %>%
  rename(Gene.names = "Gene") %>%
  left_join(hupo2, by = "Gene.names")

tovolc_protease <- filter(tovolc, Protein %in% proteases$Entry) 
write.csv(tovolc_protease, file = "results_3_humancohort_Tumor-TANM/volcano_hits_protease.csv")

tophits_protease <- slice_max(tovolc_protease, order_by = logFC, n = 10) %>%
  bind_rows(slice_min(tovolc_protease, order_by = logFC, n = 10))

simpvol2_protease <- ggplot(data = tovolc, mapping = aes(x = logFC, y = -log10(adj.P.Val), label = Gene.names)) +
  geom_point(data = filter(tovolc,  adj.P.Val > 0.01),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "black", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC > 0, adj.P.Val <= 0.01),  
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "#A52A2A", alpha = 0.4)+
  geom_point(data = filter(tovolc, logFC < 0, adj.P.Val <= 0.01), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "#4682B4", alpha = 0.4)+
  geom_text_repel(aes(label = ifelse(Protein %in% tophits_protease$Protein, as.character(Gene.names),'')),
                  hjust = 0, vjust = 0, size = 3, max.overlaps = Inf,
                  min.segment.length = 3, point.size = NA, show.legend = FALSE) +
  theme_classic() +
  #labs(title = "Tumor / TANM") +
  ylab("-log10 (adj. p-value)") +
  xlab("log2 fold change")

print(simpvol2_protease)

ggsave("results_3_humancohort_Tumor-TANM/volcano_protease.pdf", plot = simpvol2_protease, width = 63, height = 42, 
       units = "mm", dpi = 300, scale = 2)

tiff("results_3_humancohort_Tumor-TANM/volcano_protease.tiff", units = "in", width = 10, height = 6, res = 300, pointsize =(10*100/72))
print(simpvol2_protease)
dev.off()

png("results_3_humancohort_Tumor-TANM/volcano_protease.png", units = "in", width = 10, height = 6, res = 300, pointsize =(10*100/72))
print(simpvol2_protease)
dev.off()
