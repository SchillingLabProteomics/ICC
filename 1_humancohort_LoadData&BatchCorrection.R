# Batch correction of OvCa data
# v1.5 - Niko Pinter

library(sva)
library(dplyr)
library(stringr)
library(tidyr)
library(magrittr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(missMDA)
library(tibble)
library(cowplot)
library(qualpalr)
library(gridExtra)
library(randomcoloR)
library(readxl)
library(devtools)
library(diann)
library(tidyverse)

# Read and prepare DIA-NN output
dataraw <- diann_load("Data/ICC_humancohort_DIA-NNoutput.tsv")

#remove filepath from File.Name
dataraw <- mutate(dataraw, File.Name = str_replace(dataraw[[1]], ".*\\\\", ""))
dataraw <- mutate(dataraw, File.Name = str_replace(dataraw[[1]], ".raw.mzml$", ""))

#extract MQ-LFQ data from file - Genes.MQLFQ is already normalized via DIA-NN!

precursors <- diann_matrix(dataraw, pg.q = 0.01, proteotypic.only = T)
unique.genes <- diann_matrix(dataraw, id.header="Protein.Ids", 
                             quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01)
write.csv(unique.genes, "Data/ICC_humancohort_matrix.tsv")

# functions
  plot_BP <- function(input_data, title)
  {
    label_val <- round(sqrt(ncol(input_data)))
    input_data <- stack(as.data.frame(input_data))
    input_data %>%
      ggplot(aes(x = ind,
                 y = values)) +
        geom_boxplot(outlier.shape = 1) +
        ggtitle(title) +
        labs(y = "log2-transformed LFQ intensity") +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(size = 12, hjust = 0.5),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.text.x = element_text(angle = 45, size = 12, color = "black", vjust = 1, hjust = 1)) +
        scale_x_discrete(name = "", breaks = levels(input_data$ind)[seq(1, length(levels(input_data$ind)), label_val)])
  }

  plot_PCA <- function(input_data, condition, batch, title, color, show_clusters=FALSE)
  {
    if (sum(is.na(input_data)) > 0) {
      # Principal component analysis with missing values:a comparative survey of methods
      # DOI 10.1007/s11258-014-0406-z
      data <- missMDA::imputePCA(input_data)$completeObs
    } else {
      data <- input_data
    }
    pca_data <- prcomp(t(data[, seq_len(ncol(input_data))]))
    pca_data_perc <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
    df_pca_data <- data.frame(PC1 = pca_data$x[, 1],
                              PC2 = pca_data$x[, 2],
                              sample = colnames(input_data),
                              condition = condition)

    ggplot(df_pca_data, aes(PC1, PC2,
                            color = as.factor(batch),
                            shape = as.factor(condition))) +
      ggtitle(title) +
      geom_point(size = 4) +
      {
        if (show_clusters) {
          stat_ellipse(aes(PC1, PC2,
                           color = as.factor(batch),
                           group = as.factor(batch)),
                       type = "norm")
        }
      } +
      scale_color_manual(values=color) +
      theme_bw() +
      theme(legend.position = "none",
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12, color = "black")) +
      labs(x=paste0("PC1 (",pca_data_perc[1],")"),
           y=paste0("PC2 (",pca_data_perc[2],")")) +
      labs(color = "Batches", shape = "Conditions")
  }

  plot_SP <- function(input_x, input_y, samples)
  {
    batch_max <- length(samples)
    batch_random_num <- round(sqrt(batch_max))
    batch_random <- sample(samples, batch_random_num)
    sp <- list()

    for (i in seq_len(batch_random_num)) {
      sp[[i]] <- ggplot(data=data.frame(x=input_x[[batch_random[i]]],
                                        y=input_y[[batch_random[i]]]),
                        aes(x, y)) +
                  geom_point(size = 0.5) +
                  ggtitle(batch_random[i]) +
                  theme_bw() +
                  theme(legend.position = "none",
                        plot.title = element_text(size = 12, hjust = 0.5),
                        axis.title = element_blank(),
                        axis.text.y = element_text(size = 12, color = "black"),
                        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1))
    }
    gridExtra::grid.arrange(grobs = sp, ncol = 3,
                 bottom = 'input data (log2-trans. LFQ int.)',
                 left = 'batch corrected (log2-trans. LFQ int.)')
  }
############################################

# setup
###################################################################

  # set experiment/cohort title
    exp_title <- "ICC_DIA-NNoutput_batchcorrection_measurements"
    exp_title_repl <- gsub(" ", "_", exp_title, fixed = TRUE)
  ############################################

  # set max allowed NAs per protein (per batch) (in %) for ComBat -> you need to "titrate" this number -> if the number is above a certain treshold, the script crashes -> find the highest possible NA number
    expr_data_allow_NAs <- 75
  ############################################

  # calculate P.Values and adj.P.Values
    calc_p.val <- FALSE
  ############################################

  # read files
    study_design <- read_excel("Data/study_design.xlsx")
    annotation <- read_excel("Data/annotation_batch_correction.xlsx")
    expr_data <- read.csv("Data/ICC_humancohort_matrix.tsv", sep = ",")
    
    row.names(expr_data) <- expr_data$X
    expr_data <- expr_data[,-1] %>% 
      log2() %>%
      dplyr::select(-c("KB42", "KB147", "KB53", "KB213", "KB215", "KB55", "KB23", "KB209"))
    study_design <- filter(study_design, Filename %in% colnames(expr_data))
    annotation <- filter(annotation, FullRunName %in% colnames(expr_data))
    
  ############################################

  # clean study_design
    study_design <- study_design %>%
                      remove_rownames() %>%
                      column_to_rownames("Filename") %>%
                      rename(condition = "Condition") %>%
                      dplyr::select(condition)
  ############################################

  # convert annotation data to comBat format
    annotation <- annotation %>%
                    remove_rownames() %>%
                    column_to_rownames("FullRunName") %>%
                    dplyr::select(-digestion_batch) %>%
                    rename(sample = measurement_order) %>%
                    rename(batch = MS_batch) %>%
                    transform(batch = as.numeric(str_replace(batch,
                                                             "Batch_",
                                                             "")))

  ############################################

###################################################################

# merge condition from study design with annotation data
    annotation <- merge(annotation, study_design, by = 0) %>%
      column_to_rownames(var = "Row.names")
############################################

# filter rows according to expr_data_percent_NAs value
  expr_batch_n <- annotation %>%
                    group_by(batch) %>%
                    summarize(n = n()) %>%
                    mutate(max_NAs = round(n * (expr_data_allow_NAs/100)))

  for (i in seq_len(nrow(expr_batch_n))) {
   samples_per_batch <- annotation %>%
                          filter(batch == expr_batch_n$batch[[i]]) %>%
                          rownames()
   expr_batch_n$samples[i] <- data.frame(samples_per_batch)
  }

  expr_data_red_NAs <- expr_data

  k <- 0

  for (i in seq_len(nrow(expr_data))) {
    for (j in seq_len(nrow(expr_batch_n))) {
      # count NAs per batch per protein
      NAs_per_batch <- expr_data[i,] %>%
                        dplyr::select(unlist(expr_batch_n[j,]$samples)) %>%
                        is.na() %>%
                        sum()
      # drop protein if counted NAs per batch > allowed NAs per batch
      if (NAs_per_batch > expr_batch_n[j,]$max_NAs) {
        expr_data_red_NAs <- expr_data_red_NAs %>%
                              subset(!rownames(expr_data_red_NAs) %in% rownames(expr_data[i,]))
        # print(paste0("Removed ", rownames(expr_data[i,]),
        #              " (allowed NAs in batch ",j, " = ", expr_batch_n[j,]$max_NAs,
        #              " / counted NAs per batch per protein = ", NAs_per_batch))
        k <- k + 1
        # break loop if true for one of the batches per protein
        break
      }
    }
  }

  print(paste(k,"proteins removed!"))
############################################

# applying ComBat
  batch <- annotation$batch
  mod0 <- model.matrix(~1, data = annotation)
  combat_edata <- sva::ComBat(dat = expr_data_red_NAs,
                         batch = batch,
                         mod = mod0,
                         par.prior = TRUE,
                         prior.plots = FALSE)
  combat_edata <- as.data.frame(combat_edata)
############################################

# calculate P.Value and adj.P.Value
  if (calc_p.val == TRUE) {
    mod <- model.matrix(~as.factor(condition), data = annotation)

    p.val.comb <- f.pvalue(combat_edata, mod, mod0)
    p.val.adj.comb <- p.adjust(p.val.comb, method = "BH")

    p.val.comb <- data.frame(p.val.comb)
    p.val.adj.comb <- data.frame(p.val.adj.comb)

    combat_edata <- cbind(combat_edata, p.val.comb[, 1][match(rownames(combat_edata), rownames(p.val.comb))])
    colnames(combat_edata)[ncol(combat_edata)] <- "P.Value"
    combat_edata <- cbind(combat_edata, p.val.adj.comb[, 1][match(rownames(combat_edata), rownames(p.val.adj.comb))])
    colnames(combat_edata)[ncol(combat_edata)] <- "adj.P.Value"
  }
############################################

# plots
  dir.create(paste0("./results_",exp_title_repl), showWarnings = FALSE)

  color_pal <- randomcoloR::distinctColorPalette(length(levels(as.factor(annotation$batch))))

  bp_bc <- plot_BP(expr_data_red_NAs,
                   "input data")

  bp_ac <- plot_BP(combat_edata,
                   "batch corrected")

  PCA_bc <- plot_PCA(expr_data_red_NAs, annotation$condition, annotation$batch,
                     "input data", color_pal,
                     show_clusters = TRUE)

  PCA_ac <- plot_PCA(combat_edata, annotation$condition, annotation$batch,
                     "batch corrected", color_pal,
                     show_clusters = TRUE)

  sp <- plot_SP(expr_data_red_NAs, combat_edata, row.names(annotation))

  title <- cowplot::ggdraw() +
    cowplot::draw_label(exp_title, x = 0.5, y = 0.8, fontface = "bold") +
    cowplot::draw_label(paste0(expr_data_allow_NAs, "% allowed [NAs/batch]/protein"), x = 0.5, y = 0.6, size = 12) +
    cowplot::draw_label(paste0("proteins used for ComBat: ", nrow(expr_data_red_NAs), "/", nrow(expr_data)), x = 0.5, y = 0.4, size = 12) +
    cowplot::draw_label(paste0("proteins removed before ComBat: ", k, "/", nrow(expr_data)), x = 0.5, y = 0.2, size = 12)

  multi_plot_top <- cowplot::plot_grid(bp_bc, bp_ac, ncol = 2)
  multi_plot_bottom <- cowplot::plot_grid(PCA_bc, PCA_ac, ncol = 2)
  multi_plot_legend <- cowplot::get_legend(PCA_bc +
                                    guides(color = guide_legend(nrow = 2)) +
                                    theme(legend.position = "bottom",
                                          legend.title = element_text(size = 12),
                                          legend.text = element_text(size = 12)))
  multi_plot <- cowplot::plot_grid(title,
                          multi_plot_top,
                          sp,
                          multi_plot_bottom,
                          multi_plot_legend,
                          nrow = 5,
                          rel_heights = c(0.25, 1, 1, 1, 0.25))

  ggsave(paste0("./results_",exp_title_repl,"/",exp_title_repl,"_batch_corrected_max_",expr_data_allow_NAs,"%%_NAs.pdf"),
         plot = multi_plot,
         width = 30,
         height = 50,
         dpi = 72,
         units = "cm")

