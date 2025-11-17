library(tidyverse)
library(dplyr)
library(ggpubr)


# boxplots of only the candidate genes (from signature).
candidate_genes <-  c("CPNE3", "DYNC1I2", "EPS8L1", "FAXDC2", "GLS", "HECW2",   
                                "IGSF10",  "KCTD21",  "PGM2L1" ) ## approach 2: no_tmm upregulated k-fold validation.
candidate_genes <- c("HTR6", "TPGS1", "TSEN54", "WDR24", "WDR74") ## approach 2: tmm upregulated k-fold validation.


candidate_genes <- c("CALM2", "CPNE3", "DYNC1I2", "EPS8L1", "FAXDC2", "GLS", "HECW2",
                                "IGSF10", "KCTD21", "LIFR", "PGM2L1", "PYROXD1") # approach 2: no_tmm Upregulated 2-times repeat in k-fold validation.
candidate_genes <- c("DDN", "HTR6", "PIDD1", "TPGS1", 
                                 "TSEN54", "WDR24", "WDR74")  # approach 2: tmm Upregulated 2-times repeat in k-fold validation.




candidate_genes <- c("ACADM", "EPHA5", "FAXDC2", "GLIPR1L2", "HOXC9", "KIFAP3",
                                "NEK7", "PRDM2") ## approach 1: no_tmm upregulated k-fold validation.
candidate_genes <- c("CPA1", "TERT", "WDR74") ## approach 1: tmm upregulatd k-fold validation.


candidate_genes <- c("ACADM", "CPNE3", "EPHA5", "EPS8L1", "FAXDC2", "FGD4", "GLIPR1L2",
                               "HOXC9", "KIFAP3", "KYAT3", "MMP16", "NEK7",
                               "PLEKHA5", "PRDM2", "SATB1", "ZNF197") # approach 1: no_tmm Upregulated 2-times repeat in k-fold validation.
candidate_genes <- c("ALDH1A2", "CPA1",    
                                 "FOXK2", "MAGEA9", "PRR7",    
                                 "RNF126", "SPEF1", "TERT", "TSEN54", "WDR74", "XRCC3") # approach 1: tmm upregulated 2-times repeat in k-fold validation.



candidate_genes <- c("CPNE8", "EPS8L1", "FAXDC2", "IGSF10", "KIFAP3", "MYO9A",
                                "PGM2L1", "PYROXD1", "SH3GLB1") # NO_TMM upregulated: approach 3; 3-times repeat in 5-fold validation.
candidate_genes <- c("TERT", "TSEN54", "WDR24", "WDR74") # TMM upregulated: approach 3; 3-times repeat in 5-fold validation.


candidate_genes <- c("CPNE8", "EPS8L1", "FAXDC2", "IGSF10", "KIFAP3",  
                                "MYO9A", "PGM2L1", "PYROXD1", "RAB5B", "SH3GLB1") # NO_TMM upregulated: approach 3; 2 times repeat in 5-fold validation.
candidate_genes <- c("DDX39A", "TERT", "TSEN54", "WDR24", "WDR74") # TMM upregulated: approach 3; 2 times repeat in 5-fold validation.


# layout for subplots.
num_genes <- length(candidate_genes)

# Saving the plot as a PDF.
pdf("rna-seq-regression_results_noTMMup2.pdf", width = 4, height = 4)

regression_test_candidates <- lcpm[rownames(lcpm) %in% candidate_genes, ]
regression_test_candidates <- regression_test_candidates[, match(metadata_combined$SampleID, colnames(regression_test_candidates))]



for (gene in rownames(regression_test_candidates)) {
  plot_data <- data.frame(gene_Expression = as.numeric(regression_test_candidates[gene, ]),
                          TMM = metadata_combined$TMM, TMM_Case = metadata_combined$TMM_Case)
  
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMM_Case, data = plot_data)
  summary_model <- summary(model)
  r2_label <- paste0("R² = ", round(summary_model$r.squared, 3))
  
  
  # boxplot.
  p <- ggplot(plot_data, aes(x = TMM, y = gene_Expression, fill = TMM, color = TMM)) +
    geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), size = 3) +
    scale_fill_manual(values = c("ALT"="lightblue",
                                 "NO_TMM" = "lightgreen",
                                 "Telomerase"="lightpink2")) +
    scale_color_manual(values = c("ALT"="blue",
                                  "Telomerase"="darkred",
                                  "NO_TMM" = "darkgreen")) +
    theme_classic() +
    labs(title = gene, x = "Class", y = "Expression") +
    theme(
      axis.text.x = element_text(vjust = 1, hjust = 1),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    annotate("text", label = r2_label, y = max(plot_data$gene_Expression, na.rm = TRUE) * 1.05, x = 2,
             hjust = 0, size = 7) +
    stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.75) +
    stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.90)
  
  print(p)
  
  rm(summary_model)
  rm(plot_data)
  rm(model)
  rm(r2_label)
  
}

dev.off()


#######################################################################################


###########################################################################################################################

### correlation with extend score.

### correlation with expression instead of ranks.
source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")


extendScores <- RunEXTEND(as.matrix(lcpm))
telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- telomeraseScores %>%
  mutate(SampleID = ifelse(grepl("^H", SampleID),
                              gsub("\\.", "-", SampleID),
                              SampleID))

telomeraseScores <- left_join(telomeraseScores, metadata_combined[, c("SampleID", "TMM")], by = "SampleID")
telomeraseScores <- as.data.frame(telomeraseScores)
rownames(telomeraseScores) <- telomeraseScores$SampleID
telomeraseScores$SampleID = NULL


# first looking at correlation of EXTEND with TMM upregulated markers.
ranked_NBL <- apply(lcpm, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_NBL)
t_NBL <- as.data.frame(t_NBL)
t_NBL <- t_NBL[match(rownames(telomeraseScores), rownames(t_NBL)), ]

# vector to store results.
cor_values <- numeric(ncol(t_NBL))
p_values <- numeric(ncol(t_NBL))

# Looping through each gene for rank-based Spearman correlation with EXTEND.
for (i in seq_along(t_NBL)) {
  result <- cor.test(t_NBL[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# building result table for all genes for calculating FDR.
cor_results <- data.frame(
  Gene = colnames(t_NBL),
  correlation = cor_values,
  p_values = p_values
)

# calculating FDR.
cor_results$FDR <- p.adjust(cor_results$p_values, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# Now extracting only the TMM-upregulated subset.
cor_results_tmmUpregulated_extend <- cor_results[
  cor_results$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"),
]
cor_results_tmmUpregulated_extend <- cor_results_tmmUpregulated_extend %>%
  arrange(FDR)

cor_results_tmmUpregulated_extend <- cor_results[
  cor_results$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"),
]
cor_results_tmmUpregulated_extend <- cor_results_tmmUpregulated_extend %>%
  arrange(FDR)


# Now, looking at EXTEND correlated genes for NO_TMM upregulated markers.
cor_results_notmmUpregulated_extend <- cor_results[
  cor_results$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3"),
]
cor_results_notmmUpregulated_extend <- cor_results_notmmUpregulated_extend %>%
  arrange(FDR)

cor_results_notmmUpregulated_extend <- cor_results[
  cor_results$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2"),
]
cor_results_notmmUpregulated_extend <- cor_results_notmmUpregulated_extend %>%
  arrange(FDR)





######################################################################################################
######################################################################################################

## Now, instead of ranks: using expression values for TERT.

# transpose so samples are rows and genes are columns.
t_NBL <- t(lcpm)
t_NBL <- as.data.frame(t_NBL)

# vector to store results.
cor_values <- numeric(ncol(t_NBL))
p_values <- numeric(ncol(t_NBL))

# Looping through each gene for expression-based Pearson correlation with TERT.
for (i in seq_along(t_NBL)) {
  result <- cor.test(t_NBL[[i]], tert_expr, method = "pearson")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# building result table for all genes for calculating FDR.
cor_results <- data.frame(
  Gene = colnames(t_NBL),
  correlation = cor_values,
  p_values = p_values
)

# calculating FDR.
cor_results$FDR <- p.adjust(cor_results$p_values, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# Now extracting only the TMM-upregulated subset.
cor_results_tmmUpregulated_tert <- cor_results[
  cor_results$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"),
]
cor_results_tmmUpregulated_tert <- cor_results_tmmUpregulated_tert %>%
  arrange(FDR)

cor_results_tmmUpregulated_tert <- cor_results[
  cor_results$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"),
]
cor_results_tmmUpregulated_tert <- cor_results_tmmUpregulated_tert %>%
  arrange(FDR)


# Now, looking at TERT correlated genes for NO_TMM upregulated markers.
cor_results_notmmUpregulated_tert <- cor_results[
  cor_results$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3"),
]
cor_results_notmmUpregulated_tert <- cor_results_notmmUpregulated_tert %>%
  arrange(FDR)

cor_results_notmmUpregulated_tert <- cor_results[
  cor_results$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2"),
]
cor_results_notmmUpregulated_tert <- cor_results_notmmUpregulated_tert %>%
  arrange(FDR)
#############################


#######################

## correlation matrix of all genes.

library(pheatmap)
library(Hmisc)

t_NBL <- t(lcpm)
tmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", 
                                                             "MAGEA9", "SPEF1", "TERT", "WDR74"), drop = FALSE]
notmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", 
                                                               "HOXC9", "ITPRID2", "MMP16", "PRDM2") , drop = FALSE]

tmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("LCN15", "TPGS1", "TSEN54", "WDR74"), drop = FALSE]
notmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("ACADM", "CALM2", "CPNE3", "FAXDC2", 
                                                               "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") , drop = FALSE]
all_expression <- t_NBL[, colnames(t_NBL) %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", 
                                                 "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3", 
                                                 "LCN15", "TPGS1", "TSEN54", "WDR74") , drop = FALSE]

res <- rcorr(as.matrix(t_NBL), type = "pearson")

cor_matrix <- res$r     
p_matrix <- res$P


# calculating adjusted p-value.
pval_matrix_adj <- apply(p_matrix, 1, p.adjust, method = "fdr")

## now only taking the signature genes.
cor_matrix <- cor_matrix[rownames(cor_matrix) %in% colnames(all_expression), colnames(cor_matrix) %in% colnames(all_expression)]
pval_matrix_adj <- pval_matrix_adj[rownames(pval_matrix_adj) %in% colnames(all_expression), colnames(pval_matrix_adj) %in% colnames(all_expression)]


# 
# Creating a matrix of significance marks.
sig_matrix <- ifelse(pval_matrix_adj < 0.001, "***",
                     ifelse(pval_matrix_adj < 0.01, "**",
                            ifelse(pval_matrix_adj < 0.05, "*", "")))
sig_matrix[is.na(sig_matrix)] <- ""

# formatting label.
combined_labels <- matrix(
  paste0(sprintf("%.2f", cor_matrix), "\n", sig_matrix),
  nrow = nrow(cor_matrix),
  dimnames = dimnames(cor_matrix)
)

# heatmap.
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = combined_labels,  # showing significance stars + correlation value.
         number_color = "black",
         fontsize_number = 12,
         fontsize = 13,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average")

grid::grid.text("* p < 0.05, ** p < 0.01, *** p < 0.001",
                x = 0.86, y = 0.95, gp = grid::gpar(fontsize = 11))



alt_samples <- metadata_combined$SampleID[metadata_combined$TMM == "ALT"]
tel_samples <- metadata_combined$SampleID[metadata_combined$TMM == "Telomerase"]
notmm_samples <- metadata_combined$SampleID[metadata_combined$TMM == "NO_TMM"]

cor_alt <- cor((tmmUpregulated_expression[alt_samples, ]), method = "pearson")
cor_tel <- cor((tmmUpregulated_expression[tel_samples, ]), method = "pearson")
cor_noTMM <- cor((tmmUpregulated_expression[notmm_samples, ]), method = "pearson")

pheatmap(cor_tel,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         fontsize = 13,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average")



#######################################################################################

## for TARGET, looking at the correlation of telomere content vs. TMM.

t_NBL <- t(tmm_lcpm_target)

tmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("LCN15", "TPGS1", "TSEN54", "WDR74"), drop = FALSE]
notmmUpregulated_expression <- t_NBL[, colnames(t_NBL) %in%  c("ACADM", "CALM2", "CPNE3", "FAXDC2", 
                                                               "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") , drop = FALSE]
all_expression <- t_NBL[, colnames(t_NBL) %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", 
                                                 "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3", 
                                                 "LCN15", "TPGS1", "TSEN54", "WDR74") , drop = FALSE]
t_NBL <- as.data.frame(t_NBL)
t_NBL <- t_NBL[match(metadata$SampleID, rownames(t_NBL)), ]

### Now, looking at the correlation of all genes with telomere content.
telomere_content <- metadata$Telomere.Content

# vector to store results.
cor_values <- numeric(ncol(t_NBL))
p_values <- numeric(ncol(t_NBL))

# Looping through each gene for expression-based Pearson correlation with telomere content.
for (i in seq_along(t_NBL)) {
  result <- cor.test(t_NBL[[i]], telomere_content, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# building result table for all genes for calculating FDR.
cor_results <- data.frame(
  Gene = colnames(t_NBL),
  correlation = cor_values,
  p_values = p_values
)

# calculating FDR.
cor_results$FDR <- p.adjust(cor_results$p_values, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# Now extracting only the TMM-upregulated subset.
cor_results_tmmUpregulated_tc <- cor_results[
  cor_results$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"),
]
cor_results_tmmUpregulated_tc <- cor_results_tmmUpregulated_tc %>%
  arrange(FDR)

cor_results_tmmUpregulated_tc <- cor_results[
  cor_results$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"),
]
cor_results_tmmUpregulated_tc <- cor_results_tmmUpregulated_tc %>%
  arrange(FDR)


# Now, looking at NO_TMM upregulated markers.
cor_results_notmmUpregulated_tc <- cor_results[
  cor_results$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3"),
]
cor_results_notmmUpregulated_tc <- cor_results_notmmUpregulated_tc %>%
  arrange(FDR)

cor_results_notmmUpregulated_tc <- cor_results[
  cor_results$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2"),
]
cor_results_notmmUpregulated_tc <- cor_results_notmmUpregulated_tc %>%
  arrange(FDR)

# no correlation between telomere content and signatures here.
#################################################################################################

## Looking at if there is difference in telomere content between ALT, Telomerase and NO_TMM.

ggplot(metadata, aes(x = TMM, y = Telomere.Content, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "Telomere Content") +  scale_y_log10() +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


########################################################################################
#############################################################################################

# Telomere content vs. gsva per sample.
candidate_genes <- list(TMM = c("CPNE3", "DYNC1I2", "EPS8L1", "FAXDC2", "GLS", "HECW2",   
                                "IGSF10",  "KCTD21",  "PGM2L1" )) ## approach 2: no_tmm upregulated k-fold validation.
candidate_genes2 <- list(TMM = c("HTR6", "TPGS1", "TSEN54", "WDR24", "WDR74")) ## approach 2: tmm upregulated k-fold validation.


candidate_genes <- list(TMM = c("CALM2", "CPNE3", "DYNC1I2", "EPS8L1", "FAXDC2", "GLS", "HECW2",
                                "IGSF10", "KCTD21", "LIFR", "PGM2L1", "PYROXD1")) # approach 2: no_tmm Upregulated 2-times repeat in k-fold validation.
candidate_genes2 <- list(TMM = c("DDN", "HTR6", "PIDD1", "TPGS1", 
                                 "TSEN54", "WDR24", "WDR74"))  # approach 2: tmm Upregulated 2-times repeat in k-fold validation.




candidate_genes <- list(TMM = c("ACADM", "EPHA5", "FAXDC2", "GLIPR1L2", "HOXC9", "KIFAP3",
                                "NEK7", "PRDM2" )) ## approach 1: no_tmm upregulated k-fold validation.
candidate_genes2 <- list(TMM = c("CPA1", "TERT", "WDR74")) ## approach 1: tmm upregulatd k-fold validation.


candidate_genes <- list(TMM =c("ACADM", "CPNE3", "EPHA5", "EPS8L1", "FAXDC2", "FGD4", "GLIPR1L2",
                               "HOXC9", "KIFAP3", "KYAT3", "MMP16", "NEK7",
                               "PLEKHA5", "PRDM2", "SATB1", "ZNF197")) # approach 1: no_tmm Upregulated 2-times repeat in k-fold validation.
candidate_genes2 <- list(TMM = c("ALDH1A2", "CPA1",    
                                 "FOXK2", "MAGEA9", "PRR7",    
                                 "RNF126", "SPEF1", "TERT", "TSEN54", "WDR74", "XRCC3")) # approach 1: tmm upregulated 2-times repeat in k-fold validation.



candidate_genes <- list(TMM = c("CPNE8", "EPS8L1", "FAXDC2", "IGSF10", "KIFAP3", "MYO9A",
                                "PGM2L1", "PYROXD1", "SH3GLB1")) # NO_TMM upregulated: approach 3; 3-times repeat in 5-fold validation.
candidate_genes2 <- list(TMM = c("TERT", "TSEN54", "WDR24", "WDR74")) # TMM upregulated: approach 3; 3-times repeat in 5-fold validation.


candidate_genes <- list(TMM = c("CPNE8", "EPS8L1", "FAXDC2", "IGSF10", "KIFAP3",  
                                "MYO9A", "PGM2L1", "PYROXD1", "RAB5B", "SH3GLB1")) # NO_TMM upregulated: approach 3; 2 times repeat in 5-fold validation.
candidate_genes2 <- list(TMM = c("DDX39A", "TERT", "TSEN54", "WDR24", "WDR74")) # TMM upregulated: approach 3; 2 times repeat in 5-fold validation.


gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case", "Telomere.Content")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = Telomere.Content, y = GSVA_Score, colour = TMM)) +
  geom_point(size = 2, alpha = 0.7) +
 scale_x_log10() +
  labs(x = "Telomere Content",
       y = "GSVA Score") +
  theme_classic()

