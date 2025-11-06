library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(GSVA)
library(ComplexHeatmap)
library(circlize)


best_geneset <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74") # 5-fold validation tmm winner.
best_geneset2 <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") # 5-fold validation no_tmm winner.
best_geneset3 <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74",
                   "ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") # together.

best_geneset <- c("LCN15", "TPGS1", "TSEN54", "WDR74")
best_geneset2 <- c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")
best_geneset3 <- c("LCN15", "TPGS1", "TSEN54", "WDR74",
                   "ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")

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

##################################################################


# signature -- TMM upregulated.
gene_map <- bitr(best_geneset,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 8,
              readable = TRUE)
goTMMUpSignature <- go@result %>%
  filter(Count > 1)




# signature -- No_TMM upregulated.
gene_map <- bitr(best_geneset2,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

go <- groupGO(gene = gene_map$ENTREZID,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              level = 8,
              readable = TRUE)

goNoTMMUpSignature <- go@result %>%
  filter(Count > 1)


###############


# calculating correlation with known biological pathway scores (cell cycle, DNA replication, telomere maintenance).
library(msigdbr)

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

bio_pathways <- hallmark_list[c(
  # DNA repair / genome stability.
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  
  # Telomerase / proliferative signaling.
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  
  # Chromatin / plasticity
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  
  # Stress / senescence
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)]


param <- gsvaParam(as.matrix(lcpm),
                     bio_pathways,
                      kcdf = "Gaussian")

gsva_results <- gsva(param)

# for sample alignment.
signature_scores <- lcpm[rownames(lcpm) %in% best_geneset3, ]
gsva_results <- gsva_results[, colnames(signature_scores), drop = FALSE]


cor_matrix <- matrix(NA, nrow = nrow(signature_scores), ncol = nrow(gsva_results),
                     dimnames = list(rownames(signature_scores), rownames(gsva_results)))

pval_matrix <- cor_matrix

for (g in rownames(signature_scores)) {
  for (p in rownames(gsva_results)) {
    test <- cor.test(as.numeric(signature_scores[g, ]), as.numeric(gsva_results[p, ]),
                     method = "spearman")
    cor_matrix[g, p] <- test$estimate
    pval_matrix[g, p] <- test$p.value
  }
}

# Adjust p-values for multiple testing (optional, per gene or overall)
pval_matrix_adj <- apply(pval_matrix, 1, p.adjust, method = "fdr")

sig_labels <- ifelse(pval_matrix_adj < 0.001, "***",
                     ifelse(pval_matrix_adj < 0.01, "**",
                            ifelse(pval_matrix_adj < 0.05, "*", "")))

sig_labels <- t(sig_labels)

# row_order <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74",
#                   "ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")
row_order <- c("LCN15", "TPGS1", "TSEN54", "WDR74",
               "ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")
column_order <- c("HALLMARK_DNA_REPAIR", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE",
                  "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_E2F_TARGETS",
                  "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MTORC1_SIGNALING",
                  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS",
                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING",
                  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
                  "HALLMARK_INFLAMMATORY_RESPONSE")

cor_matrix <- cor_matrix[row_order, column_order]

col_categories <- c(
  "HALLMARK_DNA_REPAIR"                      = "DNA Repair",
  "HALLMARK_G2M_CHECKPOINT"                  = "DNA Repair",
  "HALLMARK_MITOTIC_SPINDLE"                 = "DNA Repair",
  "HALLMARK_P53_PATHWAY"                     = "Apoptosis",
  "HALLMARK_APOPTOSIS"                       = "Apoptosis",
  
  "HALLMARK_E2F_TARGETS"                     = "Proliferative Signaling",
  "HALLMARK_MYC_TARGETS_V1"                  = "Proliferative Signaling",
  "HALLMARK_MYC_TARGETS_V2"                  = "Proliferative Signaling",
  "HALLMARK_MTORC1_SIGNALING"                = "Proliferative Signaling",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"       = "Proliferative Signaling",
  "HALLMARK_GLYCOLYSIS"                      = "Proliferative Signaling",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Chromatin / Plasticity",
  "HALLMARK_TGF_BETA_SIGNALING"                = "Chromatin / Plasticity",
  
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "Stress / Inflammation",
  "HALLMARK_HYPOXIA"                         = "Stress / Inflammation",
  "HALLMARK_INFLAMMATORY_RESPONSE"           = "Stress / Inflammation"
)

# column split.
col_split <- factor(col_categories[column_order],
                    levels = c("DNA Repair",
                               "Apoptosis",
                               "Proliferative Signaling",
                               "Chromatin / Plasticity",
                               "Stress / Inflammation"))

# row split.
# tmm_up <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74")
# notmm_up <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")
tmm_up <- c("LCN15", "TPGS1", "TSEN54", "WDR74")
notmm_up <- c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")

row_group <- ifelse(rownames(cor_matrix) %in% tmm_up, "TMM-upregulated signature", "NO_TMM-upregulated signature")
row_split <- factor(row_group, levels = c("TMM-upregulated signature", "NO_TMM-upregulated signature"))

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))


# significance table.
sig_labels <- sig_labels[row_order, column_order, drop = FALSE]

# final ComplexHeatmap.
pdf("TMM_correlation_heatmap.pdf", width = 14, height = 10)
ht <- Heatmap(cor_matrix,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = col_split,
        row_split = row_split,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 9, fontface = "bold"),
        row_names_side = "left",
        top_annotation = HeatmapAnnotation(Category = col_split),
        left_annotation = rowAnnotation(Group = row_split,
                                        annotation_name_gp = gpar(fontsize = 9, fontface = "bold")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "black", fill = NA, lwd = 1))
          

          label <- sig_labels[i, j]
          if (label != "") {
            grid.text(label, x, y, gp = gpar(fontsize = 10, col = "black", fontface = "bold"))
          }
        },
        column_title = "Hallmark Pathways",
        row_title = "Gene Signatures",
        heatmap_legend_param = list(title = "Spearman Correlation",
                                    title_gp = gpar(fontsize = 10, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 9)))
sig_legend <- Legend(
  labels = c("*  p < 0.05", "**  p < 0.01", "***  p < 0.001"),
  type = "text",
  legend_gp = gpar(fontsize = 10),
  pch = NA,
  title = "Significance",
  title_gp = gpar(fontsize = 10, fontface = "bold")
)

draw(ht, annotation_legend_list = list(sig_legend), padding = unit(c(5, 10, 10, 10), "mm"))   
dev.off()


################################################################################################################



