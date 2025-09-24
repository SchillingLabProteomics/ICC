### Script for running multi-group limma analysis ---- 
## Miguel Cosenza 02.07.2020 & Tilman Werner

### 1. Please provide a meaningful name for the dataset/experiment ----

exper_code <- "Limma Clusters"

### 2. How many top significant hits do you want to plot (boxplots group vs intensity)? 

n_top_hits <- 20

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
library(lubridate)
library(survival)
library(survminer)
library(ggrepel)
library(readxl)
library(clusterProfiler)
`%notin%` <- Negate(`%in%`)

## Load data ----
annot_dat <- read.csv("results_4_humancohort_tumor_unsupervisedStats/annotation_cluster_k=2.csv") %>%
  dplyr::rename(Sample_ID = 'Sample') %>%
  filter(Group == "Tumor") %>%
  filter(Sample_ID %notin% c("KB106", "KB42", "KB147", "KB53", "KB213", "KB215", "KB55", "KB23", "KB209")) %>%
  arrange(Sample_ID)

expr_dat <- read.csv("results_2_humancohort_unsupervisedStats/TimsTOF_ICC_ExprMat_log2_median.csv", na = "NA") %>%
  dplyr::select(-X) %>%
  dplyr::select("Protein", annot_dat$Sample) 
expr_dat <- filter(expr_dat, rowSums(!is.na(expr_dat)) / ncol(expr_dat) >= 0.8) 

          
## Define the design matrix ----

groups <- as.factor(annot_dat$clustertest)

design <- model.matrix(~0+groups)

row.names(design) <- annot_dat$Sample_ID

colnames(design) <- colnames(design) %>% str_remove(., "groups") %>% str_trim()

# 3. DEFINE WHICH GROUPS YOU WISH TO COMPARE ----

contrast.matrix <- makeContrasts(cluster2-cluster1, levels=design) # MODIFY THIS LINE

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
            file = "results_5_humancohort_ECM-TurnoverCluster/tab_output_general_Ftest_limma.txt",
            row.names = FALSE, col.names = TRUE)

sig_hits <- dplyr::filter(output_limma2, 
                               adj.P.Val <= 0.05) %>% row.names(.)

n_significant <- length(sig_hits)

top_n_hits <- slice_min(output_limma2, n = n_top_hits, order_by = adj.P.Val, with_ties = FALSE)


#load proteome gene names

humanproteome <- read.delim("Data/uniprot-proteome UP000005640+reviewed yes(incl.GO).tab") %>%
  dplyr::rename(Protein = "Entry") 
entries       <- str_split_fixed(humanproteome$Entry.name, "_", n=Inf)
humanproteome$Entry.name <- entries[,1]


## Prep some boxplots for the top proteins with lowest P-values ----

slim_expr <- pivot_longer(expr_dat, cols = colnames(mat),
                          names_to = "Sample_ID",
                          values_to = "Abundance") 

slim_expr_g <- left_join(slim_expr, annot_dat,
                         by = "Sample_ID")

hits_expr <- filter(slim_expr_g, Protein %in% row.names(top_n_hits))
hits_expr <- left_join(hits_expr, humanproteome, by = "Protein")            


boxplots <- ggplot(hits_expr,
                   aes(x = clustertest, y = Abundance)) +
          geom_boxplot()+
          facet_wrap(Entry.name ~ .)

plot(boxplots)
## Prep volcano plots for contrasts ----

merged_limmacontr <- reshape::merge_all(list_tabular)

write.table(merged_limmacontr,
            file = "results_5_humancohort_ECM-TurnoverCluster/tab_output_per_contrasts_limma.txt",
            row.names = FALSE, col.names = TRUE)

tovolc <- merged_limmacontr %>% mutate(Differentially_expressed = case_when(adj.P.Val <= 0.05 ~ TRUE, TRUE ~ FALSE)) 
tovolc <- left_join(tovolc, humanproteome, by = "Protein")
tovolc$Gene.names <- str_split_fixed(tovolc$Gene.names, " ", n = Inf)
tovolc <- mutate(tovolc, factor = (logFC * -log10(adj.P.Val)))
write.csv(tovolc, "results_5_humancohort_ECM-TurnoverCluster/volcano_hits.csv")

tophits1 <- filter(tovolc, adj.P.Val <= 0.001) 
tophits2 <-  slice_max(tophits1, order_by = logFC, n = 15) %>%
  bind_rows(slice_min(tophits1, order_by = logFC, n = 15))

#volcanoplot 
volcanoes2 <- ggplot(data = tovolc, 
                    mapping = aes(x = logFC,
                                  y = -log10(adj.P.Val)),
                    label = Entry.name) + 
  geom_point(data = filter(tovolc, adj.P.Val > 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "black", alpha = 0.5) +
  geom_point(data = filter(tovolc, logFC < 0, adj.P.Val <= 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "blue", alpha = 0.5) +
  geom_point(data = filter(tovolc, logFC > 0, adj.P.Val <= 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "orange", alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label=ifelse(Protein %in% tophits2$Protein, as.character(Gene.names[,1]),'')),
                  hjust=0, vjust=0, size = 3, max.overlaps = Inf, min.segment.length = 2, point.size = NA) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  #labs(title = "Protein Hits", subtitle = "Differentially expressed Cluster 1 vs. Cluster 2") +
  xlab("log2 fold change") +
  ylab("-log10 adj. p-value")
volcanoes2

ggsave("results_5_humancohort_ECM-TurnoverCluster/volcano.pdf", plot = volcanoes2, width = 90, height = 73, 
       units = "mm", dpi = 300, scale = 2)

tiff("results_5_humancohort_ECM-TurnoverCluster/volcano.tiff", units = "in", width = 8, height = 6, res = 300, pointsize =(8*100/72))
volcanoes2
dev.off()


plot(volcanoes2)

### Kaplan-Meier Curve
surv_obj <- Surv(time = annot_dat$daystoevent, event = annot_dat$censored)
fit1 <- survfit(surv_obj ~ clustertest, data = annot_dat)
  
plot2 <- ggsurvplot(fit1, data = annot_dat, pval = TRUE, 
                    conf.int = FALSE, surv.median.line = "v",
                   xlab = "Days to Recurrence", ylab = "Probability", 
                   palette = c("blue", "orange"), risk.table = TRUE, pval.coord = c(1200, 0.8),
                   pval.method = FALSE)
print(plot2)

svg("results_5_humancohort_ECM-TurnoverCluster/KM.svg", width = 6, height = 5, pointsize =(5*100/72))
print(plot2)
dev.off()

plot2B <- ggsurvplot(fit1, data = annot_dat, pval = TRUE, 
                    conf.int = FALSE, 
                    xlab = "Days to Recurrence", ylab = "Probability", 
                    palette = c("blue", "orange"), risk.table = FALSE, pval.coord = c(1200, 0.8),
                    pval.method = FALSE)
print(plot2B)

svg("results_5_humancohort_ECM-TurnoverCluster/KM_graphabstr.svg", width = 4, height = 3, pointsize =(4*100/72))
print(plot2B)
dev.off()

tiff("results_5_humancohort_ECM-TurnoverCluster/KM_graphabstr.tiff", units = "in", width = 4, height = 3, res = 300, pointsize =(4*100/72))
plot2B
dev.off()

### Gene Set Enrichment
geneList <- tovolc$logFC
names(geneList) = as.character(tovolc$Protein)
geneList = sort(geneList, decreasing = TRUE)

# activate as needed
gse <- readRDS("results_5_humancohort_ECM-TurnoverCluster/GeneSetEnrichment.RData")
# gse <- gseGO(geneList, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "UNIPROT", minGSSize = 80, maxGSSize = 500,
#              pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
# 
# saveRDS(gse, file = "results_5_humancohort_ECM-TurnoverCluster/GeneSetEnrichment.RData")

ridgeplot <- ridgeplot(gse, label_format = 20) + 
  theme(axis.text.y = element_text(size = 6)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  labs(x = "log2-fold change", title = "Ontology: Biological Processes")
ridgeplot

tiff("results_5_humancohort_ECM-TurnoverCluster/ridgeplot.tiff", units = "in", width = 8, height = 6, res = 300, pointsize =(8*100/72))
ridgeplot
dev.off()

plot <- gse@result 
plot1data <- slice_max(gse@result, enrichmentScore, n = 8) %>%
  bind_rows(slice_min(gse@result, enrichmentScore, n = 8)) %>%
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
  ylab("") +
  xlab("enrichment score")
plot1

ggsave("results_5_humancohort_ECM-TurnoverCluster/dotplot.pdf", plot = plot1, width = 98, height = 54, 
       units = "mm", dpi = 300, scale = 2)

tiff("results_5_humancohort_ECM-TurnoverCluster/dotplot.tiff", units = "in", width = 8, height = 4, res = 300, pointsize =(8*100/72))
plot1
dev.off()

# EIF4A1
expr_dat2 <- read.csv("results_2_humancohort_unsupervisedStats/TimsTOF_ICC_ExprMat_log2_median.csv", na = "NA") %>%
  dplyr::select(-X) %>%
  dplyr::select("Protein", annotation$Sample) %>%
  filter(Protein == "P60842") %>%
  pivot_longer(cols = contains("KB"), names_to = "Sample", values_to = "Abundance") %>%
  left_join(annotation, by = "Sample")

plot4 <- ggplot(expr_dat2, aes(x = Tissue, y = Abundance, fill = clusterinfo)) +
  geom_line(aes(group = Patient_No)) +
  geom_boxplot(width = 0.5) +
  #geom_jitter() +
  facet_wrap(~ clustertest) +
  scale_fill_manual(values = c("cadetblue", "orange", "blue"),
                    breaks = c("TANM", "cluster2", "cluster1")) +
  stat_compare_means(comparisons = list(c("normal", "Tumor"))) +
  theme_minimal()


tiff("results_5_humancohort_ECM-TurnoverCluster/EIF4A1.tiff", units = "in", width = 4, height = 4, 
     res = 300, pointsize =(4*100/72))
plot4
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
            file = "results_5_humancohort_ECM-TurnoverCluster/KEGG_90_NM.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t"
)

setwd("results_5_humancohort_ECM-TurnoverCluster/KEGG_pathway/")
for(varhsa in hsa){
  pv.out <- pathview(gene.data =fgdata, 
                     pathway.id = varhsa, 
                     species = "hsa", 
                     out.suffix = "90perc_NM", 
                     kegg.native = TRUE,
                     same.layer = FALSE,
                     low = list(gene = "blue"), 
                     mid =list(gene = "gray"), 
                     high = list(gene = "orange"))
}


# Proteases
hupo2 <- read.delim("Data/uniprot-proteome UP000005640+reviewed yes.tab") %>%
  separate_rows(Gene.names, sep = " ")
proteases <- read_excel("Data/MEROPS_human_proteases.xlsx") %>%
  filter(!is.na(Gene)) %>%
  rename(Gene.names = "Gene") %>%
  left_join(hupo2, by = "Gene.names")

tovolc_protease <- filter(tovolc, Protein %in% proteases$Entry) 
write.csv(tovolc_protease, file = "results_5_humancohort_ECM-TurnoverCluster/volcano_hits_protease.csv")

tophits_protease <- slice_max(tovolc_protease, order_by = logFC, n = 10) %>%
  bind_rows(slice_min(tovolc_protease, order_by = logFC, n = 10)) %>%
  filter(abs(logFC) >= 0.5 & adj.P.Val <= 0.05)

simpvol2_protease <-  ggplot(data = tovolc, 
                             mapping = aes(x = logFC,
                                           y = -log10(adj.P.Val)),
                             label = Entry.name) + 
  geom_point(data = filter(tovolc, adj.P.Val > 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "black", alpha = 0.5) +
  geom_point(data = filter(tovolc, logFC < 0, adj.P.Val <= 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "blue", alpha = 0.2) +
  geom_point(data = filter(tovolc, logFC > 0, adj.P.Val <= 0.05), 
             mapping = aes(x = logFC, y = -log10(adj.P.Val)),
             color = "orange", alpha = 0.2) +
  theme_classic() +
  geom_text_repel(aes(label=ifelse(Protein %in% tophits_protease$Protein, as.character(Gene.names[,1]),'')),
                  hjust=0, vjust=0, size = 3, max.overlaps = Inf, min.segment.length = 2, point.size = NA) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  #labs(title = "Protein Hits", subtitle = "Differentially expressed Cluster 1 vs. Cluster 2") +
  xlab("log2 fold change") +
  ylab("-log10 adj. p-value")

print(simpvol2_protease)

ggsave("results_5_humancohort_ECM-TurnoverCluster/volcano_protease.pdf", plot = simpvol2_protease, width = 63, height = 42, 
       units = "mm", dpi = 300, scale = 2)

tiff("results_5_humancohort_ECM-TurnoverCluster/volcano_protease.tiff", units = "in", width = 10, height = 6, res = 300, pointsize =(10*100/72))
print(simpvol2_protease)
dev.off()

png("results_5_humancohort_ECM-TurnoverCluster/volcano_protease.png", units = "in", width = 10, height = 6, res = 300, pointsize =(10*100/72))
print(simpvol2_protease)
dev.off()
