library(EnhancedVolcano)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(biomaRt)
library(edgeR)
library(GSVA)

source("noTMM/DataCleaning.R")

# un-log the data.
gene_Expression <- round(2^Expression - 1)


# metadata <- metadata %>%
#   filter(COG.Risk.Group == "High Risk")

# filtering for protein coding genes.
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

gene_Expression <- gene_Expression[rownames(gene_Expression) %in% protein_coding_genes$hgnc_symbol, ]
gene_Expression <- gene_Expression[, colnames(gene_Expression) %in% metadata$SampleID]


# setting metadata order.
metadata <- metadata %>%
  arrange(TMM_Case, TMM)


gene_Expression <- gene_Expression[, match(metadata$SampleID, colnames(gene_Expression)), drop = FALSE]


##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

# creating differential gene expression object.
dge_TMM <- DGEList(counts=gene_Expression,group=group1)

# removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# Calculating dispersion and fitting the model.
d <- estimateDisp(dge_TMM, design, verbose=TRUE)
fit <- glmQLFit(d, design, robust = TRUE)

# contrast parameter (ALT-NO_TMM).
contrast <- makeContrasts(altVSnotmm = group1ALT - group1NO_TMM,
                          TelomeraseVSnotmm = group1Telomerase - group1NO_TMM,
                          ALTvsTelomerase = group1ALT - group1Telomerase,
                          levels = design
)


# differential expression test.
fitALT <- glmQLFTest(fit, contrast = contrast[, "altVSnotmm"])
fitTelomerase <- glmQLFTest(fit, contrast = contrast[, "TelomeraseVSnotmm"])

# results
top_ALT_target <- topTags(fitALT, n = Inf)
top_Telomerase_target <- topTags(fitTelomerase, n = Inf)



# filtering for candidate genes.
candidate_genes_ALT_target <- subset(top_ALT_target$table, FDR <= 0.05 & abs(logFC) >= 0.5)
candidate_genes_Telomerase_target <- subset(top_Telomerase_target$table, FDR <= 0.05 & abs(logFC) >= 0.5)



tmm_cpm_target  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_target <- cpm(dge_TMM, log = TRUE, prior.count = 1)   


# only including genes present in both TARGET and 0532 data.
tmm_lcpm_0532 <- read_tsv("tmm_lcpm_0532.tsv") #from 0532DataCleaning.R
tmm_lcpm_0532 <- as.data.frame(tmm_lcpm_0532)
rownames(tmm_lcpm_0532) <- tmm_lcpm_0532$...1 
tmm_lcpm_0532$...1 <- NULL

tmm_lcpm_target <- tmm_lcpm_target[rownames(tmm_lcpm_target) %in% rownames(tmm_lcpm_0532), ,drop = FALSE]
write.table(tmm_lcpm_target, file = "tmm_lcpm_target.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

candidate_genes_ALT_target <- candidate_genes_ALT_target[rownames(candidate_genes_ALT_target) %in% rownames(tmm_lcpm_target), ]
candidate_genes_Telomerase_target <- candidate_genes_Telomerase_target[rownames(candidate_genes_Telomerase_target) %in% rownames(tmm_lcpm_target), ]

#############

# t-test of the genes.
##### 2) t-test of the DEGs.

# first looking at difference between NO_TMM and Telomerase
noTMM_candidates <- tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% rownames(candidate_genes_Telomerase_target),
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("Telomerase", "NO_TMM")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "Telomerase"], drop = FALSE]
  
  
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


# Storing significant t-test results where FDR is less than 0.05.
t_test_results_Telomerase_sig <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.01, ]

#####

## t-test between NO_TMM and ALT
noTMM_candidates <- tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% rownames(candidate_genes_ALT_target),
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("ALT", "NO_TMM")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "ALT"], drop = FALSE]
  
  
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
t_test_results_ALT_sig <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.01, ]



######################################################################################################


##### 3) Linear regression test of candidate genes.

# first ALT.

regression_results_ALT <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_ALT <-  tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% t_test_results_ALT_sig$Gene,
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("ALT", "NO_TMM")],
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
regression_results_ALT_sig <- regression_results_ALT %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.3)
write.table(regression_results_ALT_sig, file = "TARGET-ALT.tsv", sep = '\t', quote = FALSE, row.names = TRUE)




### Now, turn for Telomerase.

regression_results_Telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_Telomerase <-  tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% t_test_results_Telomerase_sig$Gene,
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("Telomerase", "NO_TMM")],
  drop = FALSE]


meta_sub <- metadata[metadata$SampleID %in% colnames(dge_gene_Telomerase), ]
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

# Filtering for r-squared >= 0.2 & fdr <= 0.01.
regression_results_Telomerase_sig <- regression_results_Telomerase %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.3)
write.table(regression_results_Telomerase_sig, file = "TARGET-Telomerase.tsv", sep = '\t', quote = FALSE, row.names = TRUE)

# 
# #####################################################################################
# # boxplots of only the candidate genes.
# 
# # layout for subplots.
# num_genes <- length(candidates_tmm_upregulated_telomerase_combined$Gene)
# 
# # Saving the plot as a PDF.
# pdf("rna-seq-regression_results.pdf", width = 4, height = 4)
# 
# regression_test_candidates <- lcpm[rownames(lcpm) %in% candidates_tmm_upregulated_telomerase_combined$Gene, ]
# regression_test_candidates <- regression_test_candidates[, match(metadata_combined$SampleID, colnames(regression_test_candidates))]
# 
# 
# 
# for (gene in rownames(regression_test_candidates)) {
#   plot_data <- data.frame(gene_Expression = as.numeric(regression_test_candidates[gene, ]), 
#                           TMM = metadata_combined$TMM, TMM_Case = metadata_combined$TMM_Case)
#   
#   
#   # Fit linear model.
#   model <- lm(gene_Expression ~ TMM_Case, data = plot_data)
#   summary_model <- summary(model)
#   r2_label <- paste0("R² = ", round(summary_model$r.squared, 3))
#   
#   
#   # boxplot.
#   p <- ggplot(plot_data, aes(x = TMM, y = gene_Expression, fill = TMM, color = TMM)) +
#     geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#     geom_point(position = position_jitter(width = 0.2), size = 3) +
#     scale_fill_manual(values = c("ALT"="lightblue", 
#                                  "NO_TMM" = "lightgreen", 
#                                  "Telomerase"="lightpink2")) +
#     scale_color_manual(values = c("ALT"="blue", 
#                                   "Telomerase"="darkred", 
#                                   "NO_TMM" = "darkgreen")) +
#     theme_classic() +
#     labs(title = gene, x = "Class", y = "Expression") +
#     theme(
#       axis.text.x = element_text(vjust = 1, hjust = 1),
#       axis.title = element_text(size = 6),
#       axis.text = element_text(size = 6, face = "bold"),
#       legend.position = "none"
#     ) +
#     annotate("text", label = r2_label, y = max(plot_data$gene_Expression, na.rm = TRUE) * 1.05, x = 2,
#              hjust = 0, size = 4) +
#     stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                        method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.75) +
#     stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                        method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.90)
#   
#   print(p)
#   
#   rm(summary_model)
#   rm(plot_data)
#   rm(model)
#   rm(r2_label)
#   
# }
# 
# dev.off()

# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 

