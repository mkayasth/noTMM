library(dplyr)
library(tidyverse)
library(ggpubr)
library(edgeR)
library(ggfortify) # for pca autoplot.
library(sva)
library(caret)
library(GSVA)

source("noTMM/DataCleaning.R")
source("noTMM/0532DataCleaning.R")

# un-log the TARGET data.
gene_Expression <- round(2^Expression - 1)


common_genes <- Reduce(intersect, list(rownames(gene_Expression), rownames(geneExpression)))

counts_combined <- cbind(geneExpression[common_genes, ],
                         gene_Expression[common_genes, ])

# making metadata.
metadata_target <- metadata[, c("SampleID", "TMM", "TMM_Case", "COG.Risk.Group")]
metadata_target <- metadata_target %>%
  arrange(TMM_Case, TMM, COG.Risk.Group)

metadata_target <- metadata_target %>%
  mutate(Cohort = "TARGET")

metadata_0532_2 <- metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")]
colnames(metadata_0532_2) <- c("SampleID", "TMM", "TMM_Case")

metadata_0532_2 <- metadata_0532_2 %>%
  mutate(TMM = case_when(
    TMM == "TMM-" ~ "NO_TMM",
    TMM == "ALT+" ~ "ALT",
    TMM == "TERT+" ~ "Telomerase",
    TRUE ~ TMM))

metadata_0532_2 <- metadata_0532_2 %>%
  mutate(COG.Risk.Group = "High Risk",
         Cohort = "0532")

metadata_combined <- rbind(metadata_target, metadata_0532_2)
metadata_combined <- metadata_combined %>%
  arrange(TMM_Case, TMM, Cohort, COG.Risk.Group)

counts_combined <- counts_combined[, match(metadata_combined$SampleID, colnames(counts_combined))]

# edgeR TMM normalization.
group1 <- as.factor(metadata_combined$TMM)
design <- model.matrix(~group1+0)

dge <- DGEList(counts = counts_combined, group = group1)

dge <- calcNormFactors(dge, method = "TMM")

# cpm and logCPM normalization.
lcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Running ComBat batch correction.
batch <- metadata_combined$Cohort
mod   <- model.matrix(~ TMM, data = metadata_combined)

combat_lcpm <- ComBat(
  dat = as.matrix(lcpm),
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

### using pca to seee if batch correction worked.

# before and after ComBat.
pca_pre  <- prcomp(t(lcpm))
pca_post <- prcomp(t(combat_lcpm))

# Plot PCA colored by batch
autoplot(pca_pre,  data = metadata_combined, colour = 'Cohort') +
  ggtitle("Before ComBat (by Batch)") +
  theme_classic(base_size = 14)

autoplot(pca_post, data = metadata_combined, colour = 'Cohort') +
  ggtitle("After ComBat (by Batch)") +
  theme_classic(base_size = 14)

# PCA colored by phenotype
autoplot(pca_post, data = metadata_combined, colour = 'TMM', x = 1, y = 2) +
  ggtitle("After ComBat (by TMM)") +
  theme_classic(base_size = 14)


### combining genes significantly differentially expressed in their individual datasets to see if they still are differentially expressed in this combined dataset.
combined_ALT_target <- rbind(unique(rownames(candidate_genes_ALT_target), rownames(candidate_genes_ALT)))
combined_ALT_target <- c(combined_ALT_target)
combined_telomerase_target <- rbind(unique(rownames(candidate_genes_Telomerase_target), rownames(candidate_genes_Telomerase)))
combined_telomerase_target <- c(combined_telomerase_target)

#################################################################################

######## t-test.

# first looking at difference between NO_TMM and ALT.
noTMM_candidates <- lcpm[
  rownames(lcpm) %in% combined_ALT_target,
  colnames(lcpm) %in% metadata_combined$SampleID[metadata_combined$TMM %in% c("ALT", "NO_TMM")],
  drop = FALSE
]

# Initializing an empty data frame for storing t-test results.
t_test_results_ALT <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)



# all p-values for FDR adjustment later.
all_p_values <- numeric()

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- rownames(noTMM_candidates)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  Expression <- noTMM_candidates[rownames(noTMM_candidates) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all TMM sample expression data placed in c2_group.
  c1_group <- Expression[, colnames(Expression) %in% metadata_combined$SampleID[metadata_combined$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_combined$SampleID[metadata_combined$TMM == "ALT"], drop = FALSE]
  
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_ALT <- rbind(t_test_results_ALT, data.frame(
    Gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  
  
  # removing un-required intermediates formed while forming the t-test table above.
  rm(gene_id)
  rm(gene_Expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
}

# Calculating the FDR-adjusted p-values
t_test_results_ALT$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where FDR is less than 0.01.
t_test_results_ALT_sig_combined <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.01, ]


# Now,  looking at difference between NO_TMM and Telomerase.
noTMM_candidates <- lcpm[
  rownames(lcpm) %in% combined_telomerase_target,
  colnames(lcpm) %in% metadata_combined$SampleID[metadata_combined$TMM %in% c("Telomerase", "NO_TMM")],
  drop = FALSE
]

# Initializing an empty data frame for storing t-test results.
t_test_results_Telomerase <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)



# all p-values for FDR adjustment later.
all_p_values <- numeric()

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- rownames(noTMM_candidates)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  Expression <- noTMM_candidates[rownames(noTMM_candidates) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all TMM sample expression data placed in c2_group.
  c1_group <- Expression[, colnames(Expression) %in% metadata_combined$SampleID[metadata_combined$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_combined$SampleID[metadata_combined$TMM == "Telomerase"], drop = FALSE]
  
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_Telomerase <- rbind(t_test_results_Telomerase, data.frame(
    Gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  
  # removing un-required intermediates formed while forming the t-test table above.
  rm(gene_id)
  rm(gene_Expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
}

# Calculating the FDR-adjusted p-values
t_test_results_Telomerase$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where FDR is less than 0.01.
t_test_results_Telomerase_sig_combined <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.01, ]

##########

### Now, regression test.

##### 3) Linear regression test of candidate genes.

# first ALT.

regression_results_ALT <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_ALT <-  lcpm[
  rownames(lcpm) %in% t_test_results_ALT_sig_combined$Gene,
  colnames(lcpm) %in% metadata_combined$SampleID[metadata_combined$TMM %in% c("ALT", "NO_TMM")],
  drop = FALSE]


meta_sub <- metadata[metadata$SampleID %in% colnames(dge_gene_ALT), ]
meta_sub <- meta_sub[match(colnames(dge_gene_ALT), meta_sub$SampleID), ]


# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene_ALT)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene_ALT[gene, ]), 
                     TMMstatus = meta_sub$TMM)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_ALT <- rbind(regression_results_ALT, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # removing intermediates.
  rm(data)
  rm(summary_model)
  rm(gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# adjusted p-values.
regression_results_ALT$Adj.P.Value <- p.adjust(regression_results_ALT$p_value, method = "fdr")

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_ALT_sig_combined <- regression_results_ALT %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.3)
write.table(regression_results_ALT_sig_combined, file = "Combined-ALT.tsv", sep = '\t', quote = FALSE, row.names = TRUE)


##### Now, turn for Telomerase.

regression_results_Telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_Telomerase <-  lcpm[
  rownames(lcpm) %in% t_test_results_Telomerase_sig_combined$Gene,
  colnames(lcpm) %in% metadata_combined$SampleID[metadata_combined$TMM %in% c("Telomerase", "NO_TMM")],
  drop = FALSE]


meta_sub <- metadata_combined[metadata_combined$SampleID %in% colnames(dge_gene_Telomerase), ]
meta_sub <- meta_sub[match(colnames(dge_gene_Telomerase), meta_sub$SampleID), ]


# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene_Telomerase)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene_Telomerase[gene, ]), 
                     TMMstatus = meta_sub$TMM)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_Telomerase <- rbind(regression_results_Telomerase, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # removing intermediates.
  rm(data)
  rm(summary_model)
  rm(gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# adjusted p-values.
regression_results_Telomerase$Adj.P.Value <- p.adjust(regression_results_Telomerase$p_value, method = "fdr")

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_Telomerase_sig_combined <- regression_results_Telomerase %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.3)
write.table(regression_results_Telomerase_sig_combined, file = "combined-Telomerase.tsv", sep = '\t', quote = FALSE, row.names = TRUE)

############## 

candidates_tmm_upregulated_telomerase_combined <- regression_results_Telomerase_sig_combined %>%
  filter(estimate > 0)

candidates_notmm_upregulated_telomerase_combined <- regression_results_Telomerase_sig_combined %>%
  filter(estimate < 0)

candidates_tmm_upregulated_ALT_combined <- regression_results_ALT_sig_combined %>%
  filter(estimate < 0)

candidates_notmm_upregulated_ALT_combined <- regression_results_ALT_sig_combined %>%
  filter(estimate > 0)



candidates_tmm_upregulated_ALT_combined <- candidates_tmm_upregulated_ALT_combined %>%
  arrange(desc(R.Squared))
candidates_tmm_upregulated_telomerase_combined <- candidates_tmm_upregulated_telomerase_combined %>%
  arrange(desc(R.Squared))
candidates_notmm_upregulated_ALT_combined <- candidates_notmm_upregulated_ALT_combined %>%
  arrange(desc(R.Squared))
candidates_notmm_upregulated_telomerase_combined <- candidates_notmm_upregulated_telomerase_combined %>%
  arrange(desc(R.Squared))


### finding the best signature.

##### dividing lcpm into training and testing set.
set.seed(123)
lcpm <- lcpm[, match(metadata_combined$SampleID, colnames(lcpm))]
metadata_combined$Strata <- interaction(metadata_combined$TMM, metadata_combined$Cohort, drop = TRUE)
train_idx <- createDataPartition(metadata_combined$Cohort, p = 0.8, list = FALSE)

lcpm_train <- lcpm[, train_idx, drop = FALSE]
lcpm_test <- lcpm[, -train_idx, drop = FALSE]

metadata_train <- metadata_combined[train_idx, ]
metadata_test <- metadata_combined[-train_idx, ]


candidates_tmm_upregulated_combined <- rbind(candidates_tmm_upregulated_telomerase_combined, candidates_tmm_upregulated_ALT_combined)
a <- intersect(candidates_tmm_upregulated_telomerase_combined$Gene, candidates_tmm_upregulated_ALT_combined$Gene)



tmm_target_upregulated <- run_exhaustive_forward_auc(expr_matrix = lcpm_train, metadata = metadata_train, 
                                                     candidate_genes = unique(candidates_tmm_upregulated_combined$Gene), 
                                                     phenotype_col = "TMM_Case",
                                                     label_one = "TMM", label_two = "NO_TMM",
                                                     max_genes = 20,
                                                     pivot_gene = "WDR74")



candidates_notmm_upregulated_combined <- rbind(candidates_notmm_upregulated_ALT_combined, candidates_notmm_upregulated_telomerase_combined)
b <- intersect(candidates_notmm_upregulated_telomerase_combined$Gene, candidates_notmm_upregulated_ALT_combined$Gene)

tmm_target_downregulated <- run_exhaustive_forward_auc(expr_matrix = lcpm_train, metadata = metadata_train, 
                                                       candidate_genes = unique(candidates_notmm_upregulated_combined$Gene), 
                                                       phenotype_col = "TMM_Case",
                                                       label_one = "NO_TMM", label_two = "TMM",
                                                       max_genes = 20,
                                                       pivot_gene = "FAXDC2")



##### Looking at the boxplots for the signature genes.

# First, on training set.

candidate_genes <- list(TMM = c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5",
                                "CRYBG2", "NBPF6", "RUVBL1", "GALR2", "DDX39A", "C6orf132", "C1QTNF4")) #0.8 with wdr74 as pivot.

# candidate_genes <- list(TMM = c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5", "CRYBG2", "NBPF6", "C1QTNF4", "DPEP3")) #0.7 with wdr74 as pivot.

# candidate_genes <- list(TMM = c("WDR74", "LCN15", "ALG1L2", "TERT", "CPA1", "EEF1D", "OTX2")) # 0.6 with wdr74 as pivot.

gsvapar <- gsvaParam(as.matrix(lcpm_train), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_train[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


# Now, on testing set.
gsvapar <- gsvaParam(as.matrix(lcpm_test), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_test[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

## Now, trying on target data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## Now, trying on 0532 data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# Now, trying on Ackerman data.
gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat
ackerman_NB <- log2(ackerman_NB + 1)

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
# ackerman_metadata <- ackerman_metadata %>%
#   filter(Risk == "YES")

ackerman_metadata <- ackerman_metadata %>%
  filter(
    !(TERTRearrangement == "+" & TMM_Category != "Telomerase"),
    
    !( (ATRXMutation != "-" & ATRXMutation != "<NA>" & !is.na(ATRXMutation))
       & TMM_Category != "ALT"),
    
    !(APB == "+" & TMM_Category != "ALT")
  )

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case, TMM_Category)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

#############################

##### Looking at the boxplots for the signature genes -- now, downregulated genes.


# candidate_genes2 <- list(TMM = c("FAXDC2", "STRADB", "KLHL2", "PRDM2", "AGL", "OSBPL8",
#                                  "ITPRID2", "PAK1", "EPHA5", "FAM162B")) # 0.7 equally divided FAXDC2 pivot.

# candidate_genes2 <- list(TMM = c("FAXDC2", "MYO5A", "STRADB", "KLHL2", "MOB1B", "PRDM2", "EXTL2",
#                                  "MAPK1", "FAM162B", "PPP3CB")) # 0.6 equally divided FAXDC2 pivot.

# 
candidate_genes2 <- list(TMM = c("FAXDC2", "PDP1", "ACADM", "ITPRID2", "DDAH1", "SLC10A7", "DYNC1I2", "FGD4", "SATB1",
                                 "KLHL2", "HOXC9", "STRADB", "EIF4G3", "OSBPL8", "KCTD21", "EPHA5", "ETV2", "GNPTAB")) # 0.8 equally divided FAXDC2 pivot.

gsvapar <- gsvaParam(as.matrix(lcpm_train), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_train[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


# Now, on testing set.
gsvapar <- gsvaParam(as.matrix(lcpm_test), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_test[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

## Now, trying on target data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## Now, trying on 0532 data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# Now, trying on Ackerman data.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


##### subtracting TMM upregulated genes by downregulated genes and making boxplots.


## trying it in training set.
gsvapar <- gsvaParam(as.matrix(lcpm_train), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(lcpm_train), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_train[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

## doing for test set.
gsvapar <- gsvaParam(as.matrix(lcpm_test), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(lcpm_test), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_test[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## doing for target.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

### doing the same for 0532 data.

gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)




#### doing the same thing for Ackerman.
gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat
ackerman_NB <- log2(ackerman_NB + 1)

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
# ackerman_metadata <- ackerman_metadata %>%
#   filter(Risk == "YES")

ackerman_metadata <- ackerman_metadata %>%
  filter(
    !(TERTRearrangement == "+" & TMM_Category != "Telomerase"),
    
    !( (ATRXMutation != "-" & ATRXMutation != "<NA>" & !is.na(ATRXMutation))
       & TMM_Category != "ALT"),
    
    !(APB == "+" & TMM_Category != "ALT")
  )

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case, TMM_Category)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

# genes positively upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "gsva_Score_down"
  )

# genes positively upregulated in NO_TMM.
gsvapar2 <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL


gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "gsva_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$gsva_Score_up - gsva_long$gsva_Score_down

ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


######
# k-fold validation.

tmm_upregulated_kfold <- run_kfold_auc(expr_matrix = lcpm, metadata = metadata_combined, candidate_genes = candidates_tmm_upregulated_combined$Gene,
                                       phenotype_col = "TMM_Case", label_one = "TMM", label_two = "NO_TMM",
                                       max_genes = 20, k = 5,
                                       pivot_gene = "WDR74")

notmm_upregulated_kfold <- run_kfold_auc(expr_matrix = lcpm, metadata = metadata_combined, candidate_genes = candidates_notmm_upregulated_combined$Gene,
                                         phenotype_col = "TMM_Case", label_one = "NO_TMM", label_two = "TMM",
                                         max_genes = 20, k = 5,
                                         pivot_gene = "FAXDC2")

## genes present in at least three of the five folds -- for tmm upregulated genes.
all_genes <- unlist(lapply(tmm_upregulated_kfold$fold_results, function(x) x$selected_genes))
gene_counts <- table(all_genes)

kfold_tmm_upregulated_genes <- names(gene_counts[gene_counts >= 3])


## genes present in at least three of the five folds --  for no_tmm upregulated genes.
all_genes <- unlist(lapply(notmm_upregulated_kfold$fold_results, function(x) x$selected_genes))
gene_counts <- table(all_genes)

kfold_notmm_upregulated_genes <- names(gene_counts[gene_counts >= 3])

candidate_genes <- list(TMM = c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74"))
candidate_genes2 <- list(TMM = c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2"))

# upregulated in NO_TMM.
gsvapar <- gsvaParam(as.matrix(lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_combined[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## Now, trying on target data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## Now, trying on 0532 data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# Now, trying on Ackerman data.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

#############################

##### Looking at the boxplots for the signature genes -- now, downregulated genes.

gsvapar <- gsvaParam(as.matrix(lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_combined[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

## Now, trying on target data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## Now, trying on 0532 data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# Now, trying on Ackerman data.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## incorporating both GSVAs.
gsvapar <- gsvaParam(as.matrix(lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_combined[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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


## doing for target.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

### doing the same for 0532 data.

gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2",
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred",
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

#### doing the same thing for Ackerman.
gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat
ackerman_NB <- log2(ackerman_NB + 1)

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
# ackerman_metadata <- ackerman_metadata %>%
#   filter(Risk == "YES")

ackerman_metadata <- ackerman_metadata %>%
  filter(
    !(TERTRearrangement == "+" & TMM_Category != "Telomerase"),
    
    !( (ATRXMutation != "-" & ATRXMutation != "<NA>" & !is.na(ATRXMutation))
       & TMM_Category != "ALT"),
    
    !(APB == "+" & TMM_Category != "ALT")
  )

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case, TMM_Category)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

# genes positively upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "gsva_Score_down"
  )

# genes positively upregulated in NO_TMM.
gsvapar2 <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL


gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "gsva_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$gsva_Score_up - gsva_long$gsva_Score_down

ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2",
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred",
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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



###############################################################################################



