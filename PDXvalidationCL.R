library(dplyr)
library(tidyverse)
library(biomaRt)
library(edgeR)
library(GSVA)
library(ggpubr)

# validating results in cell Line.
cellLineExpression <- readRDS("Kallisto_CellLines_BC_Counts_Final_06222023.RDS")
cellLineMetadata <- read.delim("CLs_Metadata_ForBatchCorrection_11042025.txt", sep = '\t')
cellLineMetadata <- cellLineMetadata[cellLineMetadata$TMM %in% c("ALT", "TERT", "EST"), ]
cellLineMetadata <- cellLineMetadata %>%
  arrange(TMM)

colnames(cellLineExpression) <- gsub("\\.", "-", colnames(cellLineExpression))

# filtering for protein coding genes.
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

cellLineExpression <- cellLineExpression[cellLineExpression$GeneSymbol %in% protein_coding_genes$hgnc_symbol, ]

rownames(cellLineExpression) <- cellLineExpression$GeneSymbol
cellLineExpression <- cellLineExpression[, colnames(cellLineExpression) %in% cellLineMetadata$RNAseq_SampleID]
cellLineExpression <- cellLineExpression[, match(cellLineMetadata$RNAseq_SampleID, colnames(cellLineExpression))]

##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(cellLineMetadata$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

# creating differential gene expression object.
dge_TMM <- DGEList(counts=cellLineExpression,group=group1)

# removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")

tmm_cpm_cl  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_cl <- cpm(dge_TMM, log = TRUE, prior.count = 1)


## boxplot to see whats up.
best_geneset <- list(TMM = c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74")) # 5-fold validation tmm winner.
best_geneset2 <- list(TMM = c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")) # 5-fold validation no_tmm winner.

best_geneset <- list(TMM = c("LCN15", "TPGS1", "TSEN54", "WDR74")) # 5-fold validation tmm winner.
best_geneset2 <- list(TMM = c("ACADM", "CALM2", "CPNE3", "FAXDC2",              # 5-fold validation no_tmm winner.
                             "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3"))

gsvapar <- gsvaParam(as.matrix(tmm_lcpm_cl), best_geneset, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )


gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_cl), best_geneset2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, cellLineMetadata[, c("RNAseq_SampleID", "TMM")], by = c("SampleID" = "RNAseq_SampleID"))
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT" = "lightpink2",
                               "EST" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("TERT"="darkred",
                                "EST" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT","EST")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT","EST")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


#####################################################################################################################

# validating results in PDX.

pdxExpression <- readRDS("KallistoGE_BC_Counts_Final_06202023.RDS")
pdxMetadata <- read.delim("PDXs_Metadata_11042025.txt", sep = '\t')
pdxMetadata <- pdxMetadata[pdxMetadata$TMM %in% c("ALT", "TERT", "EST"), ]
pdxMetadata <- pdxMetadata %>%
  arrange(TMM)

colnames(pdxExpression) <- gsub("\\.", "-", colnames(pdxExpression))

# filtering for protein coding genes.
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

pdxExpression <- pdxExpression[pdxExpression$GeneSymbol %in% protein_coding_genes$hgnc_symbol, ]

rownames(pdxExpression) <- pdxExpression$GeneSymbol
pdxExpression <- pdxExpression[, colnames(pdxExpression) %in% pdxMetadata$RNAseq_SampleID]
pdxExpression <- pdxExpression[, match(pdxMetadata$RNAseq_SampleID, colnames(pdxExpression))]

##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(pdxMetadata$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

# creating differential gene expression object.
dge_TMM <- DGEList(counts=pdxExpression,group=group1)

# removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")

tmm_cpm_pdx  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_pdx <- cpm(dge_TMM, log = TRUE, prior.count = 1)


## boxplot to see whats up.
best_geneset <- list(TMM = c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74")) # 5-fold validation tmm winner.
best_geneset2 <- list(TMM = c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")) # 5-fold validation no_tmm winner.

best_geneset <- list(TMM = c("LCN15", "TPGS1", "TSEN54", "WDR74")) # 5-fold validation tmm winner.
best_geneset2 <- list(TMM = c("ACADM", "CALM2", "CPNE3", "FAXDC2",              # 5-fold validation no_tmm winner.
                              "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3"))

gsvapar <- gsvaParam(as.matrix(tmm_lcpm_pdx), best_geneset, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )


gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_pdx), best_geneset2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, pdxMetadata[, c("RNAseq_SampleID", "TMM")], by = c("SampleID" = "RNAseq_SampleID"))
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT" = "lightpink2",
                               "EST" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("TERT"="darkred",
                                "EST" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT","EST")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("TERT","EST")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

