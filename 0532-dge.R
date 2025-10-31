#############


source("noTMM/0532DataCleaning.R")

# t-test of the genes.
##### 2) t-test of the DEGs.

# first looking at difference between NO_TMM and Telomerase
noTMM_candidates <- tmm_lcpm[
  rownames(tmm_lcpm) %in% rownames(candidate_genes_Telomerase),
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("TERT+", "TMM-")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TMM-"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TERT+"], drop = FALSE]
  
  
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
t_test_results_Telomerase_sig <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.01, ]

#####

## t-test between NO_TMM and ALT
noTMM_candidates <- tmm_lcpm[
  rownames(tmm_lcpm) %in% rownames(candidate_genes_ALT),
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("ALT+", "TMM-")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TMM-"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "ALT+"], drop = FALSE]
  
  
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
t_test_results_ALT_sig <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.05, ]



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

dge_gene_ALT <-  tmm_lcpm[
  rownames(tmm_lcpm) %in% t_test_results_ALT_sig$Gene,
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("ALT+", "TMM-")],
  drop = FALSE]


meta_sub <- metadata_0532[metadata_0532$RNAseq_SampleID %in% colnames(dge_gene_ALT), ]
meta_sub <- meta_sub[match(colnames(dge_gene_ALT), meta_sub$RNAseq_SampleID), ]


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
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.2)

regression_results_ALT_sig <- regression_results_ALT_sig %>%
  arrange(desc(R.Squared))




### Now, turn for Telomerase.

regression_results_Telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_Telomerase <-  tmm_lcpm[
  rownames(tmm_lcpm) %in% t_test_results_Telomerase_sig$Gene,
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("TERT+", "TMM-")],
  drop = FALSE]


meta_sub <- metadata_0532[metadata_0532$RNAseq_SampleID %in% colnames(dge_gene_Telomerase), ]
meta_sub <- meta_sub[match(colnames(dge_gene_Telomerase), meta_sub$RNAseq_SampleID), ]


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
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.2)

regression_results_Telomerase_sig <- regression_results_Telomerase_sig %>%
  arrange(desc(R.Squared))




