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
  # DNA repair pathways.
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_G2M_CHECKPOINT",
 
  # apoptosis pathways 
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  
  # Telomerase / proliferative signaling.
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  
  # Cell motility signaling.
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_APICAL_SURFACE",
  #"REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  #"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  #"REACTOME_RHO_GTPASE_CYCLE",
  #"REACTOME_ANTIGEN_PRESENTATION",
  
  
  # oxidative and metabolic stress response.
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  
  # inflammatory and cytokine stress.
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  
)]


param <- gsvaParam(as.matrix(lcpm),
                     bio_pathways,
                      kcdf = "Gaussian")

gsva_results <- gsva(param)

# doing correlation for all the genes so that adj-p value makes sense.
pathway_scores <- gsva_results
pathway_scores <- pathway_scores[, colnames(lcpm), drop = FALSE]

# initializing result matrices.
cor_matrix <- matrix(NA, nrow = nrow(lcpm), ncol = nrow(pathway_scores),
                     dimnames = list(rownames(lcpm), rownames(pathway_scores)))
pval_matrix <- cor_matrix

# computing correlation.
for (p in rownames(pathway_scores)) {
  cor_tests <- apply(lcpm, 1, function(x)
    cor.test(as.numeric(x), as.numeric(pathway_scores[p, ]), method = "spearman"))
  
  cor_matrix[, p] <- sapply(cor_tests, `[[`, "estimate")
  pval_matrix[, p] <- sapply(cor_tests, `[[`, "p.value")
}

# adjusting p-values for multiple testing per pathway.
pval_matrix_adj <- apply(pval_matrix, 2, p.adjust, method = "fdr")

# now just filtering our signature genes.
cor_matrix <- cor_matrix[rownames(cor_matrix) %in% best_geneset3, ]
pval_matrix_adj <- pval_matrix_adj[rownames(pval_matrix_adj) %in% best_geneset3, ]

sig_labels <- ifelse(pval_matrix_adj < 0.001, "***",
                     ifelse(pval_matrix_adj < 0.01, "**",
                            ifelse(pval_matrix_adj < 0.05, "*", "")))



row_order <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74",
                  "ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")
row_order <- c("LCN15", "TPGS1", "TSEN54", "WDR74",
               "ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")
column_order <- c("HALLMARK_DNA_REPAIR", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_E2F_TARGETS", 
                  "HALLMARK_MITOTIC_SPINDLE","HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", 
                  "HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_APICAL_JUNCTION",
  "HALLMARK_APICAL_SURFACE","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE")

cor_matrix <- cor_matrix[row_order, column_order]

col_categories <- c(
  "HALLMARK_DNA_REPAIR" = "DNA repair pathways",
  "HALLMARK_G2M_CHECKPOINT" = "DNA repair pathways",
  "HALLMARK_P53_PATHWAY" = "Apoptosis",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_E2F_TARGETS" = "Proliferative signaling",
  "HALLMARK_MITOTIC_SPINDLE" = "Proliferative signaling",
  "HALLMARK_MYC_TARGETS_V1" = "Proliferative signaling",
  "HALLMARK_MYC_TARGETS_V2"= "Proliferative signaling",
  "HALLMARK_MTORC1_SIGNALING"= "Proliferative signaling",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"= "Proliferative signaling",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"= "Proliferative signaling",
  "HALLMARK_GLYCOLYSIS" = "Proliferative signaling",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Cell Motility & Plasticity",
  "HALLMARK_TGF_BETA_SIGNALING" = "Cell Motility & Plasticity",
  "HALLMARK_APICAL_JUNCTION" = "Cell Motility & Plasticity",
  "HALLMARK_APICAL_SURFACE" = "Cell Motility & Plasticity",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "Oxidative & Metabolic Stress Response",
  "HALLMARK_HYPOXIA" = "Oxidative & Metabolic Stress Response",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Oxidative & Metabolic Stress Response",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "Inflammatory & Cytokine Response",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "Inflammatory & Cytokine Response",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory & Cytokine Response",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "Inflammatory & Cytokine Response"
)

# column split.
col_split <- factor(col_categories[column_order],
                    levels = c("DNA repair pathways",
                               "Apoptosis",
                               "Proliferative signaling",
                               "Cell Motility & Plasticity",
                               "Oxidative & Metabolic Stress Response",
                               "Inflammatory & Cytokine Response"))

category_colors <- c(
  "DNA repair pathways" = "blue",       
  "Apoptosis" = "green",                          
  "Proliferative signaling" = "red",            
  "Cell Motility & Plasticity" = "yellow",           
  "Oxidative & Metabolic Stress Response" = "orange",
  "Inflammatory & Cytokine Response" = "brown"
)

# row split.
tmm_up <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74")
notmm_up <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2")
tmm_up <- c("LCN15", "TPGS1", "TSEN54", "WDR74")
notmm_up <- c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")

row_group <- ifelse(rownames(cor_matrix) %in% tmm_up, "TMM-upregulated signature", "NO_TMM-upregulated signature")
row_split <- factor(row_group, levels = c("TMM-upregulated signature", "NO_TMM-upregulated signature"))

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))


# significance table.
sig_labels <- sig_labels[row_order, column_order, drop = FALSE]

# final ComplexHeatmap.
pdf("TMM_correlation_heatmap2.pdf", width = 14, height = 10)
ht <- Heatmap(cor_matrix,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = col_split,
        row_split = row_split,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 9, fontface = "bold"),
        row_names_side = "left",
        top_annotation = HeatmapAnnotation(Category = col_split, col = list(Category = category_colors)),
        left_annotation = rowAnnotation(" " = row_split,
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

draw(ht, annotation_legend_list = list(sig_legend), merge_legends = TRUE, padding = unit(c(5, 10, 10, 10), "mm"))   
dev.off()


################################################################################################################



