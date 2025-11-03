##### clusterProfiler.

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)


best_geneset <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74") # 5-fold validation tmm winner.
best_geneset2 <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") # 5-fold validation no_tmm winner.

# noTMM upregulated gene.
gene_map <- bitr(candidates_notmm_upregulated_combined$Gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 5,
              readable = TRUE)
goNoTMMUp <- go@result %>%
  filter(Count > 0)


# TMM upregulated gene.
gene_map <- bitr(candidates_tmm_upregulated_combined$Gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 5,
              readable = TRUE)
goTMMUp <- go@result %>%
  filter(Count > 0)


# signature -- TMM upregulated.
gene_map <- bitr(best_geneset,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 5,
              readable = TRUE)
goTMMUpSignature <- go@result %>%
  filter(Count > 0)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 6,
              readable = TRUE)
goTMMUpSignature <- go@result %>%
  filter(Count > 0)



# signature -- No_TMM upregulated.
gene_map <- bitr(best_geneset2,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 5,
              readable = TRUE)

goNoTMMUpSignature <- go@result %>%
  filter(Count > 0)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 6,
              readable = TRUE)
goNoTMMUpSignature <- go@result %>%
  filter(Count > 0)







###############




# calculating correlation with known biological pathway scores (cell cycle, DNA replication, telomere maintenance).
library(msigdbr)

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

bio_pathways <- hallmark_list[c("HALLMARK_G2M_CHECKPOINT",
                                "HALLMARK_E2F_TARGETS",
                                "HALLMARK_DNA_REPAIR",
                                "HALLMARK_TELOMERE_MAINTENANCE",
                                "HALLMARK_MITOTIC_SPINDLE")]

param <- gsvaParam(as.matrix(lcpm),
                     list(bio_pathways),
                      kcdf = "Gaussian")




gsva_results <- gsva(param, verbose = TRUE)

# quantify functional relatedness between the genes based on GO semantic similarity.
library(GOSemSim)


# Convert gene symbols to Entrez IDs
gene_ids <- bitr(best_geneset2, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# Compute semantic similarity matrix for Biological Process (BP)
sim_matrix <- mgeneSim(gene_ids, semData = godata("org.Hs.eg.db", ont = "BP"), measure = "Wang")

gene_map <- bitr(gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Replace row/column names
rownames(sim_matrix) <- gene_map$SYMBOL[match(rownames(sim_matrix), gene_map$ENTREZID)]
colnames(sim_matrix) <- gene_map$SYMBOL[match(colnames(sim_matrix), gene_map$ENTREZID)]

# Cluster and visualize
heatmap(sim_matrix)
mean(sim_matrix[upper.tri(sim_matrix)])


