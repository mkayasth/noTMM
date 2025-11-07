library(tidyverse)
library(dplyr)
library(ggpubr)


# boxplots of only the candidate genes (from signature).

candidate_genes <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74")
candidate_genes <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")

candidate_genes <- c("LCN15", "TPGS1", "TSEN54", "WDR74") #up in TMM.
candidate_genes <- c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") #up in no_tmm.

# layout for subplots.
num_genes <- length(candidate_genes)

# Saving the plot as a PDF.
pdf("rna-seq-regression_results_TMMup.pdf", width = 4, height = 4)

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

# Looking at associations with TERT/EXTEND.
source("noTMM/combinedDataset.R")

candidates_tmm_upregulated_combined 
candidates_notmm_upregulated_combined 
 
# first looking at correlation of TERT ranks with TMM upregulated markers.
ranked_NBL <- apply(lcpm, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_NBL)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT

tmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_tmm_upregulated_combined$Gene, drop = FALSE]
tmmUpregulated_ranks <- tmmUpregulated_ranks[match(rownames(t_NBL), rownames(tmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(tmmUpregulated_ranks))
p_values <- numeric(ncol(tmmUpregulated_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(tmmUpregulated_ranks)) {
  result <- cor.test(tmmUpregulated_ranks[[i]], tert_expr, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_tmmUpregulated <- data.frame(
  Gene = colnames(tmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_tmmUpregulated$FDR <- p.adjust(cor_results_tmmUpregulated$pvalue, method = "fdr")
cor_results_tmmUpregulated <- cor_results_tmmUpregulated %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_tmmUpregulated <- cor_results_tmmUpregulated %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

#############################

# Now, looking at correlation of TERT ranks with NO_TMM upregulated markers.
ranked_NBL <- apply(lcpm, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_NBL)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT

notmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_notmm_upregulated_combined$Gene, drop = FALSE]
notmmUpregulated_ranks <- notmmUpregulated_ranks[match(rownames(t_NBL), rownames(notmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(notmmUpregulated_ranks))
p_values <- numeric(ncol(notmmUpregulated_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(notmmUpregulated_ranks)) {
  result <- cor.test(notmmUpregulated_ranks[[i]], tert_expr, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_notmmUpregulated <- data.frame(
  Gene = colnames(notmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_notmmUpregulated$FDR <- p.adjust(cor_results_notmmUpregulated$pvalue, method = "fdr")
cor_results_notmmUpregulated <- cor_results_notmmUpregulated %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_notmmUpregulated <- cor_results_notmmUpregulated %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

###########################################################################################################################

cor_results_notmmUpregulated[cor_results_notmmUpregulated$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") , ]
cor_results_tmmUpregulated[cor_results_tmmUpregulated$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"), ]

cor_results_notmmUpregulated[cor_results_notmmUpregulated$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", 
                                                                      "HECW2", "IGSF10", "KIF13A", "KIFAP3") , ]
cor_results_tmmUpregulated[cor_results_tmmUpregulated$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74") , ]

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

tmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_tmm_upregulated_combined$Gene, drop = FALSE]
tmmUpregulated_ranks <- tmmUpregulated_ranks[match(rownames(telomeraseScores), rownames(tmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(tmmUpregulated_ranks))
p_values <- numeric(ncol(tmmUpregulated_ranks))

# Looping through each gene for rank-based Spearman correlation with EXTEND.
for (i in seq_along(tmmUpregulated_ranks)) {
  result <- cor.test(tmmUpregulated_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_tmmUpregulated_extend <- data.frame(
  Gene = colnames(tmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_tmmUpregulated_extend$FDR <- p.adjust(cor_results_tmmUpregulated_extend$pvalue, method = "fdr")
cor_results_tmmUpregulated_extend <- cor_results_tmmUpregulated_extend %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_tmmUpregulated <- cor_results_tmmUpregulated_extend %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

##### Now, looking at EXTEND correlated genes for NO_TMM upregulated markers.
notmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_notmm_upregulated_combined$Gene, drop = FALSE]
notmmUpregulated_ranks <- notmmUpregulated_ranks[match(rownames(telomeraseScores), rownames(notmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(notmmUpregulated_ranks))
p_values <- numeric(ncol(notmmUpregulated_ranks))

# Looping through each gene for rank-based Spearman correlation with EXTEND.
for (i in seq_along(notmmUpregulated_ranks)) {
  result <- cor.test(notmmUpregulated_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_notmmUpregulated_extend <- data.frame(
  Gene = colnames(notmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_notmmUpregulated_extend$FDR <- p.adjust(cor_results_notmmUpregulated_extend$pvalue, method = "fdr")
cor_results_notmmUpregulated_extend <- cor_results_notmmUpregulated_extend %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_notmmUpregulated <- cor_results_notmmUpregulated_extend %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))


#################################

cor_results_notmmUpregulated_extend[cor_results_notmmUpregulated_extend$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") , ]
cor_results_tmmUpregulated_extend[cor_results_tmmUpregulated_extend$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"), ]

cor_results_notmmUpregulated_extend[cor_results_notmmUpregulated_extend$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") , ]
cor_results_tmmUpregulated_extend[cor_results_tmmUpregulated_extend$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"), ]

##################################

######################################################################################################
######################################################################################################

## Now, instead of ranks: using expression values.

# transpose so samples are rows and genes are columns.
t_NBL <- t(lcpm)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT

tmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_tmm_upregulated_combined$Gene, drop = FALSE]
tmmUpregulated_ranks <- tmmUpregulated_ranks[match(rownames(t_NBL), rownames(tmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(tmmUpregulated_ranks))
p_values <- numeric(ncol(tmmUpregulated_ranks))

# Looping through each gene for expression-based Spearman correlation with TERT.
for (i in seq_along(tmmUpregulated_ranks)) {
  result <- cor.test(tmmUpregulated_ranks[[i]], tert_expr, method = "pearson")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_tmmUpregulated <- data.frame(
  Gene = colnames(tmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_tmmUpregulated$FDR <- p.adjust(cor_results_tmmUpregulated$pvalue, method = "fdr")
cor_results_tmmUpregulated <- cor_results_tmmUpregulated %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_tmmUpregulated <- cor_results_tmmUpregulated %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

#############################

# Now, looking at correlation of TERT expression with NO_TMM upregulated markers.

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_NBL)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT

notmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_notmm_upregulated_combined$Gene, drop = FALSE]
notmmUpregulated_ranks <- notmmUpregulated_ranks[match(rownames(t_NBL), rownames(notmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(notmmUpregulated_ranks))
p_values <- numeric(ncol(notmmUpregulated_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(notmmUpregulated_ranks)) {
  result <- cor.test(notmmUpregulated_ranks[[i]], tert_expr, method = "pearson")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_notmmUpregulated <- data.frame(
  Gene = colnames(notmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_notmmUpregulated$FDR <- p.adjust(cor_results_notmmUpregulated$pvalue, method = "fdr")
cor_results_notmmUpregulated <- cor_results_notmmUpregulated %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_notmmUpregulated <- cor_results_notmmUpregulated %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

###########################

cor_results_notmmUpregulated[cor_results_notmmUpregulated$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") , ]
cor_results_tmmUpregulated[cor_results_tmmUpregulated$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"), ]

cor_results_notmmUpregulated[cor_results_notmmUpregulated$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") , ]
cor_results_tmmUpregulated[cor_results_tmmUpregulated$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"), ]

##########################


### correlation of expression values with extend score.

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

# transpose so samples are rows and genes are columns.
t_NBL <- t(lcpm)
t_NBL <- as.data.frame(t_NBL)

tmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_tmm_upregulated_combined$Gene, drop = FALSE]
tmmUpregulated_ranks <- tmmUpregulated_ranks[match(rownames(telomeraseScores), rownames(tmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(tmmUpregulated_ranks))
p_values <- numeric(ncol(tmmUpregulated_ranks))

# Looping through each gene for expression-based Pearson correlation with EXTEND.
for (i in seq_along(tmmUpregulated_ranks)) {
  result <- cor.test(tmmUpregulated_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "pearson")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_tmmUpregulated_extend <- data.frame(
  Gene = colnames(tmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_tmmUpregulated_extend$FDR <- p.adjust(cor_results_tmmUpregulated_extend$pvalue, method = "fdr")
cor_results_tmmUpregulated_extend <- cor_results_tmmUpregulated_extend %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_tmmUpregulated <- cor_results_tmmUpregulated_extend %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

##### Now, looking at EXTEND correlated genes for NO_TMM upregulated markers.
notmmUpregulated_ranks <- t_NBL[, colnames(t_NBL) %in% candidates_notmm_upregulated_combined$Gene, drop = FALSE]
notmmUpregulated_ranks <- notmmUpregulated_ranks[match(rownames(telomeraseScores), rownames(notmmUpregulated_ranks)), ]

# vector to store results.
cor_values <- numeric(ncol(notmmUpregulated_ranks))
p_values <- numeric(ncol(notmmUpregulated_ranks))

# Looping through each gene for expression-based Pearson correlation with EXTEND.
for (i in seq_along(notmmUpregulated_ranks)) {
  result <- cor.test(notmmUpregulated_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "pearson")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results_notmmUpregulated_extend <- data.frame(
  Gene = colnames(notmmUpregulated_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results_notmmUpregulated_extend$FDR <- p.adjust(cor_results_notmmUpregulated_extend$pvalue, method = "fdr")
cor_results_notmmUpregulated_extend <- cor_results_notmmUpregulated_extend %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_notmmUpregulated <- cor_results_notmmUpregulated_extend %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))


##########################

cor_results_notmmUpregulated_extend[cor_results_notmmUpregulated_extend$Gene %in% c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") , ]
cor_results_tmmUpregulated_extend[cor_results_tmmUpregulated_extend$Gene %in% c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"), ]

cor_results_notmmUpregulated[cor_results_notmmUpregulated$Gene %in% c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3") , ]
cor_results_tmmUpregulated[cor_results_tmmUpregulated$Gene %in% c("LCN15", "TPGS1", "TSEN54", "WDR74"), ]
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







