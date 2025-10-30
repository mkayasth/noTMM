library(tidyverse)
library(dplyr)
library(ggpubr)
library(biomaRt)
library(cmapR)
library(GSVA)
library(edgeR)
library(uwot)

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


colnames(geneExpression) <- gsub("\\.", "-", colnames(geneExpression))
geneExpression <- geneExpression[, match(metadata_0532$RNAseq_SampleID, colnames(geneExpression)), drop = FALSE]



##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata_0532$TMM)

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
contrast <- makeContrasts(altVSnotmm = group1ALT. - group1TMM.,
                          TelomeraseVSnotmm = group1TERT. - group1TMM.,
                          ALTvsTelomerase = group1ALT. - group1TERT.,
                          levels = design
)


# differential expression test.
fitALT <- glmQLFTest(fit, contrast = contrast[, "altVSnotmm"])
fitTelomerase <- glmQLFTest(fit, contrast = contrast[, "TelomeraseVSnotmm"])

# results
top_ALT <- topTags(fitALT, n = Inf)
top_Telomerase <- topTags(fitTelomerase, n = Inf)



# filtering for candidate genes.
candidate_genes_ALT <- subset(top_ALT$table, round(FDR, 2) <= 0.05 & abs(logFC) >= 0.5)
candidate_genes_Telomerase <- subset(top_Telomerase$table, round(FDR, 2) <= 0.05 & abs(logFC) >= 0.5)



tmm_cpm  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
write.table(tmm_lcpm, file = "tmm_lcpm_0532.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

################################################################


## doing this after making DGE-targetRawCounts.R.
# only including genes present in TARGET raw counts.
tmm_lcpm_target <- read_tsv("tmm_lcpm_target.tsv") # from DGE-targetRawCounts.R ::)
tmm_lcpm_target <- as.data.frame(tmm_lcpm_target)
rownames(tmm_lcpm_target) <- tmm_lcpm_target$...1 
tmm_lcpm_target$...1 <- NULL

tmm_lcpm <- tmm_lcpm[rownames(tmm_lcpm) %in% rownames(tmm_lcpm_target), ,drop = FALSE]
candidate_genes_ALT <- candidate_genes_ALT[rownames(candidate_genes_ALT) %in% rownames(tmm_lcpm), ]
candidate_genes_Telomerase <- candidate_genes_Telomerase[rownames(candidate_genes_Telomerase) %in% rownames(tmm_lcpm), ]

###################################################################################














