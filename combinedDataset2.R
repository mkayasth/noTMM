### Combining two datasets and two candidate gene set.

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



##### doing NO_TMM vs. TMM in TARGET and NO_TMM vs. TMM in 0532 and seeing what comes up.

###############################################

library(EnhancedVolcano)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(biomaRt)
library(edgeR)
library(GSVA)

### First, working with TARGET set.

source("noTMM/DataCleaning.R")

# un-log the data.
gene_Expression <- round(2^Expression - 1)


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
group1 <- as.factor(metadata$TMM_Case)

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
contrast <- makeContrasts(notmmVStmm = group1NO_TMM - group1TMM,
                          levels = design
)

# differential expression test.
fitTMM <- glmQLFTest(fit, contrast = contrast[, "notmmVStmm"])

# results
top_notmm_target <- topTags(fitTMM, n = Inf)


# filtering for candidate genes.
candidate_genes_tmm_target <- subset(top_notmm_target$table, FDR <= 0.05 & abs(logFC) >= 0.5)

tmm_cpm_target  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_target <- cpm(dge_TMM, log = TRUE, prior.count = 1)   


# only including genes present in both TARGET and 0532 data.
tmm_lcpm_0532 <- read_tsv("tmm_lcpm_0532.tsv") #from 0532DataCleaning.R
tmm_lcpm_0532 <- as.data.frame(tmm_lcpm_0532)
rownames(tmm_lcpm_0532) <- tmm_lcpm_0532$...1 
tmm_lcpm_0532$...1 <- NULL

tmm_lcpm_target <- tmm_lcpm_target[rownames(tmm_lcpm_target) %in% rownames(tmm_lcpm_0532), ,drop = FALSE]


candidate_genes_tmm_target <- candidate_genes_tmm_target[rownames(candidate_genes_tmm_target) %in% rownames(tmm_lcpm_target), ]

#############

# t-test of the genes.
##### 2) t-test of the DEGs.

noTMM_candidates <- tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% rownames(candidate_genes_tmm_target),
 ,
  drop = FALSE
]

# Initializing an empty data frame for storing t-test results.
t_test_results_target <- data.frame(
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
  c1_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM_Case == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM_Case == "TMM"], drop = FALSE]
  
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_target <- rbind(t_test_results_target, data.frame(
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
t_test_results_target$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where FDR is less than 0.05.
t_test_results_target_sig <- t_test_results_target[round(t_test_results_target$fdr_t_test, 2) <= 0.01, ]



######################################################################################################


##### 3) Linear regression test of candidate genes.

regression_results_target <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_target <-  tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% t_test_results_target_sig$Gene,
  ,
  drop = FALSE]




# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene_target)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene_target[gene, ]), 
                     TMMstatus = metadata$TMM_Case)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_target <- rbind(regression_results_target, data.frame(
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
regression_results_target$Adj.P.Value <- p.adjust(regression_results_target$p_value, method = "fdr")

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_target_sig <- regression_results_target %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.3)

##########################################################################################################
#################################################################################################################


##### Now, running the same for 0532 data.

geneExpression <- readRDS("Kallisto_PTs_BC_Counts_Final_06222023.RDS")
metadata_0532 <- read_delim("TMM_052_MetaFinal_09162025.txt")



# filtering for protein coding genes.
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

geneExpression <- geneExpression[geneExpression$GeneSymbol %in% protein_coding_genes$hgnc_symbol, ]

# changing GeneSymbol column to rownames.
rownames(geneExpression) <- geneExpression$GeneSymbol
geneExpression$GeneSymbol = NULL

# setting metadata order.
metadata_0532 <- metadata_0532 %>%
  arrange(TMMCase, TMM)

# changing metadata_0532 colnames to match other data.
metadata_0532 <- metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")]
colnames(metadata_0532) <- c("SampleID", "TMM", "TMM_Case")
metadata_0532 <- metadata_0532 %>%
  mutate(TMM = case_when(
    TMM == "TMM-" ~ "NO_TMM",
    TMM == "ALT+" ~ "ALT",
    TMM == "TERT+" ~ "Telomerase",
    TRUE ~ TMM))


colnames(geneExpression) <- gsub("\\.", "-", colnames(geneExpression))
geneExpression <- geneExpression[, match(metadata_0532$SampleID, colnames(geneExpression)), drop = FALSE]


##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata_0532$TMM_Case)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

colnames(design) <- make.names(colnames(design))
colnames(design)


# creating differential gene expression object.
dge_TMM <- DGEList(counts=geneExpression,group=group1)

# removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# Calculating dispersion and fitting the model.
d <- estimateDisp(dge_TMM, design, verbose=TRUE)
fit <- glmQLFit(d, design, robust = TRUE)

# contrast parameter (ALT-NO_TMM).
contrast <- makeContrasts(tmmVSnotmm = group1NO_TMM - group1TMM,
                          levels = design
)

# differential expression test.
fit0532 <- glmQLFTest(fit, contrast = contrast[, "tmmVSnotmm"])

# results
top_0532 <- topTags(fit0532, n = Inf)

# filtering for candidate genes.
candidate_genes_0532 <- subset(top_0532$table, round(FDR, 2) <= 0.05 & abs(logFC) >= 0.5)

tmm_cpm  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

#################################################################################################




########################################################################################
