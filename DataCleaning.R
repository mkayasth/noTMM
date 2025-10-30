###
# raw counts data from TARGET NB xenahub.
# fpkm not very appropriate for DGE because of normalization of sequence depth and gene length within sample, but factor not comparable between samples.
# proceeding with log(raw counts + 1) data.
###


library(limma)
library(tidyverse)
library(biomaRt)
library(dplyr)

# Loading bulk-seq log counts RNA and metadata file.
Expression <- read_delim("TARGET-NBL.star_counts.tsv")
ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
metadata <- read_delim("Metadata_TARGETFinal_08012024.txt")

# removing Ensembl ID data after . in the name in the counts data.
Expression$Ensembl_ID <- gsub("\\..*", "", Expression$Ensembl_ID)

# filtering for protein coding genes (not done).
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

# Expression <- Expression %>%
#   filter(Ensembl_ID %in% protein_coding_genes$ensembl_gene_id)

### columns that contain expression data -- we will replace Ensembl ID with hgnc symbol as rowname now.
expr_data <- Expression[, setdiff(colnames(Expression), "Ensembl_ID")]

# Calculating proportion of zeros per row and filtering rows/genes with >=90% zeros.
zero_fraction <- rowSums(expr_data == 0) / ncol(expr_data)
Expression <- Expression[zero_fraction < 0.90, ]

# Annotating Ensembl -> GeneSymbol.
annotation <- getBM(filters = 'ensembl_gene_id',
                    attributes= c("ensembl_gene_id",
                                  "hgnc_symbol"),
                    values = Expression$Ensembl_ID,
                    mart = mart)


# removing duplicates after mapping ensembl_id to hgnc_symbol so that hgnc_symbol can be set as rownames.
annotation <- annotation[!duplicated(annotation$hgnc_symbol) & !duplicated(annotation$hgnc_symbol, fromLast = TRUE), ]

Expression <- Expression[Expression$Ensembl_ID %in% annotation$ensembl_gene_id, ]
annotation <- annotation[match(Expression$Ensembl_ID, annotation$ensembl_gene_id), ]


# Setting rownames to HGNC symbols.
Expression <- as.data.frame(Expression)
rownames(Expression) <- annotation$hgnc_symbol
Expression$Ensembl_ID = NULL

## only selecting samples in metadata and vice-versa. Also, only selecting genes also present in fpkm set.
Expression <- Expression[rownames(Expression) %in% rownames(ExpressionFpkm), ]

colnames(Expression) <- gsub("-", ".", colnames(Expression))
Expression <- Expression[, colnames(Expression) %in% metadata$SampleID]
metadata <- as.data.frame(metadata[metadata$SampleID %in% colnames(Expression), ])



###
metadata_ALT <- metadata[metadata$TMM %in% c("ALT", "NO_TMM"), ]
metadata_Telomerase <- metadata[metadata$TMM %in% c("Telomerase", "NO_TMM"), ]


metadata_ALT <- metadata_ALT %>%
  arrange(TMM_Case, TMM)
ExpressionALT <- Expression[, colnames(Expression) %in% metadata_ALT$SampleID]
ExpressionALT <- ExpressionALT[, match(metadata_ALT$SampleID, colnames(ExpressionALT))]

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM_Case, TMM)
ExpressionTelomerase <- Expression[, colnames(Expression) %in% metadata_Telomerase$SampleID]
ExpressionTelomerase <- ExpressionTelomerase[, match(metadata_Telomerase$SampleID, colnames(ExpressionTelomerase))]

