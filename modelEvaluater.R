library(GSVA)
library(pROC)
library(caret)
library(tidyverse)
library(dplyr)

################################################################################################
###############################################################################################################


evaluate_gsva_classifier <- function(expr_train, expr_test,
                                     metadata_train, metadata_test,
                                     geneset,
                                     phenotype_col = "TMM_Case",
                                     label_one = "TMM", label_two = "NO_TMM"
) {
  
  
  # gsva computation.
  combined <- cbind(expr_train, expr_test)
  valid_genes <- intersect(geneset, rownames(combined))
  if (length(valid_genes) == 0)
    stop("None of the genes in the gene set are in the expression matrices.")
  
  param <- gsvaParam(exprData = combined,
                     geneSets = list(sig = valid_genes),
                     kcdf = "Gaussian")
  mat <- gsva(param, verbose = FALSE)
  gsva_vec <- as.numeric(mat["sig", ])
  gsva_ids <- colnames(mat)
  
  # scaling gsva scores.
  gsva_vec <- as.numeric(scale(gsva_vec))
  
  # combining with metadata.
  metadata_combined <- rbind(metadata_train, metadata_test)
  metadata_aligned <- metadata_combined[match(gsva_ids, metadata_combined$SampleID), ]
  stopifnot(all(gsva_ids == metadata_aligned$SampleID))
  
  gsva_df <- data.frame(SampleID = gsva_ids,
                        GSVA_score = gsva_vec,
                        Phenotype = metadata_aligned[[phenotype_col]])
  
  # splitting GSVA scores.
  train_data <- subset(gsva_df, SampleID %in% metadata_train$SampleID)
  test_data  <- subset(gsva_df, SampleID %in% metadata_test$SampleID)
  message("train_data:", nrow(train_data), "  test_data:", nrow(test_data))
  
  # ROC training.
  roc_train <- roc(response = train_data$Phenotype,
                   predictor = train_data$GSVA_score,
                   levels = c(label_two, label_one),
                   direction = "<")
  
  coords_df <- as.data.frame(coords(roc_train, "best", best.method = "youden",
                                    ret = c("threshold", "sensitivity", "specificity")))
  threshold <- as.numeric(coords_df$threshold[1])
  if (is.na(threshold) || is.infinite(threshold)) {
    warning("Threshold returned Inf/NA; using median of training scores.")
    threshold <- median(train_data$GSVA_score, na.rm = TRUE)
  }
  message("Chosen threshold: ", round(threshold, 5))
  
  # test evaluation.
  pred <- ifelse(test_data$GSVA_score > threshold, label_one, label_two)
  stopifnot(length(pred) == nrow(test_data))
  
  cm <- confusionMatrix(factor(pred, levels = c(label_one, label_two)),
                        factor(test_data$Phenotype, levels = c(label_one, label_two)))
  
  roc_test <- roc(response = test_data$Phenotype,
                  predictor = test_data$GSVA_score,
                  levels = c(label_two, label_one),
                  direction = "<")
  
  # results.
  results <- list(
    train_auc = as.numeric(auc(roc_train)),
    test_auc  = as.numeric(auc(roc_test)),
    threshold = threshold,
    metrics = data.frame(
      Accuracy  = cm$overall["Accuracy"],
      Precision = cm$byClass["Precision"],
      Recall    = cm$byClass["Recall"],
      F1        = cm$byClass["F1"]
    ),
    confusion = cm$table,
    gsva_scores = gsva_df
  )
  return(results)
}

set.seed(123)
lcpm <- lcpm[, match(metadata_combined$SampleID, colnames(lcpm))]
metadata_combined$Strata <- interaction(metadata_combined$TMM, metadata_combined$Cohort, drop = TRUE)
train_idx <- createDataPartition(metadata_combined$Cohort, p = 0.8, list = FALSE)

lcpm_train <- lcpm[, train_idx, drop = FALSE]
lcpm_test <- lcpm[, -train_idx, drop = FALSE]

metadata_train <- metadata_combined[train_idx, ]
metadata_test <- metadata_combined[-train_idx, ]


best_geneset <- c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5",
                  "CRYBG2", "NBPF6", "RUVBL1", "GALR2", "DDX39A", "C6orf132", "C1QTNF4") #0.8 split.
best_geneset <- c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5", "CRYBG2", "NBPF6", "C1QTNF4", "DPEP3") # 0.7 split.
best_geneset <- c("WDR74", "LCN15", "ALG1L2", "TERT", "CPA1", "EEF1D", "OTX2") # 0.6 split.


best_geneset <- c("FAXDC2", "PDP1", "ACADM", "ITPRID2", "DDAH1", "SLC10A7", "DYNC1I2", "FGD4",
                  "SATB1", "KLHL2", "HOXC9", "STRADB", "EIF4G3", "OSBPL8") # 0.8 split.
best_geneset <- c("FAXDC2", "STRADB", "KLHL2", "PRDM2", "AGL", "OSBPL8", "ITPRID2", "PAK1", "EPHA5", "FAM162B") # 0.7 split.
best_geneset <- c("FAXDC2", "MYO5A", "STRADB", "KLHL2", "MOB1B", "PRDM2", "EXTL2", "MAPK1", "FAM162B", "PPP3CB") # 0.6 split.

best_geneset <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74") # 5-fold validation tmm winner.
best_geneset <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") # 5-fold validation no_tmm winner.

res <- evaluate_gsva_classifier(
  expr_train = lcpm_train,
  expr_test = lcpm_test,
  metadata_train = metadata_train,
  metadata_test = metadata_test,
  geneset = best_geneset,
  phenotype_col = "TMM_Case",
  label_one = "TMM", label_two = "NO_TMM"
)

# Inspect results
res$train_auc
res$test_auc
res$metrics
res$confusion


######################################################################################

evaluate_gsva_external <- function(expr_new, metadata_new,
                                   geneset,
                                   phenotype_col = "TMM_Case",
                                   label_one = "TMM", label_two = "NO_TMM",
                                   threshold,
                                   kcdf = "Gaussian") {
  
  expr_new <- as.matrix(expr_new)
  
  # computing gsva.
  valid_genes <- intersect(geneset, rownames(expr_new))
  if (length(valid_genes) == 0)
    stop("None of the genes in the gene set are in the expression matrix.")
  
  param <- gsvaParam(exprData = expr_new,
                     geneSets = list(sig = valid_genes),
                     kcdf = kcdf)
  mat <- gsva(param, verbose = FALSE)
  gsva_vec <- as.numeric(mat["sig", ])
  gsva_ids <- colnames(mat)
  
  # scaling gsva scores.
  gsva_vec <- as.numeric(scale(gsva_vec))
  
  metadata_aligned <- metadata_new[match(gsva_ids, metadata_new$SampleID), ]
  stopifnot(all(gsva_ids == metadata_aligned$SampleID))
  
  gsva_df <- data.frame(SampleID = gsva_ids,
                        GSVA_score = gsva_vec,
                        Phenotype = metadata_aligned[[phenotype_col]])
  
  # applying threshold for classification.
  pred <- ifelse(gsva_df$GSVA_score > threshold, label_one, label_two)
  
  # confusion matrices.
  cm <- confusionMatrix(factor(pred, levels = c(label_one, label_two)),
                        factor(gsva_df$Phenotype, levels = c(label_one, label_two)))
  
  # 
  roc_ext <- roc(response = gsva_df$Phenotype,
                 predictor = gsva_df$GSVA_score,
                 levels = c(label_two, label_one),
                 direction = "<")
  
  # result.
  results <- list(
    auc = as.numeric(auc(roc_ext)),
    threshold_used = threshold,
    metrics = data.frame(
      Accuracy  = cm$overall["Accuracy"],
      Precision = cm$byClass["Precision"],
      Recall    = cm$byClass["Recall"],
      F1        = cm$byClass["F1"]
    ),
    confusion = cm$table,
    gsva_scores = gsva_df
  )
  
  return(results)
}

###
source("noTMM/DataCleaning.R")
source("noTMM/0532DataCleaning.R")

metadata_0532 <- metadata_0532[, c("TMM", "TMMCase", "RNAseq_SampleID")]
colnames(metadata_0532) <- c("TMM", "TMM_Case", "SampleID")

gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')

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

ackerman_NB_log2 <- log2(ackerman_NB + 1)

###
res2 <- evaluate_gsva_external(expr_new = tmm_lcpm_target, metadata_new = metadata,
                               geneset = best_geneset,
                               phenotype_col = "TMM_Case",
                               label_one = "NO_TMM", label_two = "TMM",
                               threshold =  0.62802,
                               kcdf = "Gaussian")
res2



###########################################################################################


##############################


## difference between gsva of upregulated and downregulated gene sets.

evaluate_gsva_classifier <- function(expr_train, expr_test,
                                     metadata_train, metadata_test,
                                     geneset_tmmUp, geneset_tmmDown,
                                     phenotype_col = "TMM_Case",
                                     label_one = "NO_TMM", label_two = "TMM"
) {
  
  
  # gsva computation.
  combined <- cbind(expr_train, expr_test)
  valid_tmmUp <- intersect(geneset_tmmUp, rownames(combined))
  valid_tmmDown <- intersect(geneset_tmmDown, rownames(combined))
  
  if (length(valid_tmmUp) == 0 && length(valid_tmmDown) == 0)
    stop("None of the genes in the gene set are in the expression matrices.")
  
  # computing gsva score for both gene sets.
  
  param <- gsvaParam(exprData = combined,
                     geneSets = list(TMM_Up = valid_tmmUp,
                                     NO_TMM_Up = valid_tmmDown),
                     kcdf = "Gaussian")
  mat <- gsva(param, verbose = FALSE)
  
  
  gsva_vec <- as.numeric(mat["NO_TMM_Up", ] - mat["TMM_Up", ])
  gsva_ids <- colnames(mat)
  
  # scaling gsva scores.
  gsva_vec <- as.numeric(scale(gsva_vec))
  
  # combining with metadata.
  metadata_combined <- rbind(metadata_train, metadata_test)
  metadata_aligned <- metadata_combined[match(gsva_ids, metadata_combined$SampleID), ]
  stopifnot(all(gsva_ids == metadata_aligned$SampleID))
  
  gsva_df <- data.frame(SampleID = gsva_ids,
                        GSVA_score = gsva_vec,
                        Phenotype = metadata_aligned[[phenotype_col]])
  
  # splitting GSVA scores.
  train_data <- subset(gsva_df, SampleID %in% metadata_train$SampleID)
  test_data  <- subset(gsva_df, SampleID %in% metadata_test$SampleID)
  message("train_data:", nrow(train_data), "  test_data:", nrow(test_data))
  
  # ROC training.
  roc_train <- roc(response = train_data$Phenotype,
                   predictor = train_data$GSVA_score,
                   levels = c(label_two, label_one),
                   direction = "<")
  
  coords_df <- as.data.frame(coords(roc_train, "best", best.method = "youden",
                                    ret = c("threshold", "sensitivity", "specificity")))
  threshold <- as.numeric(coords_df$threshold[1])
  if (is.na(threshold) || is.infinite(threshold)) {
    warning("Threshold returned Inf/NA; using median of training scores.")
    threshold <- median(train_data$GSVA_score, na.rm = TRUE)
  }
  message("Chosen threshold: ", round(threshold, 5))
  
  # test evaluation.
  pred <- ifelse(test_data$GSVA_score > threshold, label_one, label_two)
  stopifnot(length(pred) == nrow(test_data))
  
  cm <- confusionMatrix(factor(pred, levels = c(label_one, label_two)),
                        factor(test_data$Phenotype, levels = c(label_one, label_two)))
  
  roc_test <- roc(response = test_data$Phenotype,
                  predictor = test_data$GSVA_score,
                  levels = c(label_two, label_one),
                  direction = "<")
  
  # results.
  results <- list(
    train_auc = as.numeric(auc(roc_train)),
    test_auc  = as.numeric(auc(roc_test)),
    threshold = threshold,
    metrics = data.frame(
      Accuracy  = cm$overall["Accuracy"],
      Precision = cm$byClass["Precision"],
      Recall    = cm$byClass["Recall"],
      F1        = cm$byClass["F1"]
    ),
    confusion = cm$table,
    gsva_scores = gsva_df
  )
  return(results)
}

set.seed(123)
lcpm <- lcpm[, match(metadata_combined$SampleID, colnames(lcpm))]
metadata_combined$Strata <- interaction(metadata_combined$TMM, metadata_combined$Cohort, drop = TRUE)
train_idx <- createDataPartition(metadata_combined$Cohort, p = 0.8, list = FALSE)

lcpm_train <- lcpm[, train_idx, drop = FALSE]
lcpm_test <- lcpm[, -train_idx, drop = FALSE]

metadata_train <- metadata_combined[train_idx, ]
metadata_test <- metadata_combined[-train_idx, ]


best_geneset1 <- c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5",
                   "CRYBG2", "NBPF6", "RUVBL1", "GALR2", "DDX39A", "C6orf132", "C1QTNF4") #0.8 split.
best_geneset1 <- c("WDR74", "USH1G", "TERT", "ALG1L2", "SLC1A5", "CRYBG2", "NBPF6", "C1QTNF4", "DPEP3") # 0.7 split.
best_geneset <- c("WDR74", "LCN15", "ALG1L2", "TERT", "CPA1", "EEF1D", "OTX2") # 0.6 split.


best_geneset2 <- c("FAXDC2", "PDP1", "ACADM", "ITPRID2", "DDAH1", "SLC10A7", "DYNC1I2", "FGD4",
                   "SATB1", "KLHL2", "HOXC9", "STRADB", "EIF4G3", "OSBPL8") # 0.8 split.
best_geneset2 <- c("FAXDC2", "STRADB", "KLHL2", "PRDM2", "AGL", "OSBPL8", "ITPRID2", "PAK1", "EPHA5", "FAM162B") # 0.7 split.
best_geneset <- c("FAXDC2", "MYO5A", "STRADB", "KLHL2", "MOB1B", "PRDM2", "EXTL2", "MAPK1", "FAM162B", "PPP3CB") # 0.6 split.

best_geneset1 <- c("ALG1L2", "ALOX12B", "CPA1", "DDX39A", "MAGEA9", "SPEF1", "TERT", "WDR74") # 5-fold validation tmm winner.
best_geneset2 <- c("ACADM", "EIF4G3", "EPS8L1", "FAXDC2", "FGD4", "HOXC9", "ITPRID2", "MMP16", "PRDM2") # 5-fold validation no_tmm winner.

best_geneset1 <- c("LCN15", "TPGS1", "TSEN54", "WDR74")
best_geneset2 <- c("ACADM", "CALM2", "CPNE3", "FAXDC2", "GLS", "HECW2", "IGSF10", "KIF13A", "KIFAP3")

res <- evaluate_gsva_classifier(expr_train = lcpm_train, 
                                expr_test = lcpm_test,
                                metadata_train = metadata_train, metadata_test = metadata_test,
                                geneset_tmmUp = best_geneset1, geneset_tmmDown = best_geneset2,
                                phenotype_col = "TMM_Case",
                                label_one = "NO_TMM", label_two = "TMM"
)



evaluate_gsva_external <- function(expr_new, metadata_new,
                                   geneset_tmmUp, geneset_tmmDown,
                                   phenotype_col = "TMM_Case",
                                   label_one = "TMM", label_two = "NO_TMM",
                                   threshold,
                                   kcdf = "Gaussian") {
  
  expr_new <- as.matrix(expr_new)
  
  # computing gsva.
  valid_tmmUp <- intersect(geneset_tmmUp, rownames(expr_new))
  valid_tmmDown <- intersect(geneset_tmmDown, rownames(expr_new))
  
  if (length(valid_tmmUp) == 0 && length(valid_tmmDown) == 0)
    stop("None of the genes in the gene set are in the expression matrices.")
  
  param <- gsvaParam(exprData = expr_new,
                     geneSets = list(TMM_Up = valid_tmmUp,
                                     NO_TMM_Up = valid_tmmDown),
                     kcdf = "Gaussian")
  mat <- gsva(param, verbose = FALSE)
  gsva_vec <- as.numeric(mat["NO_TMM_Up", ] - mat["TMM_Up", ])
  gsva_ids <- colnames(mat)
  
  # scaling gsva scores.
  gsva_vec <- as.numeric(scale(gsva_vec))
  
  metadata_aligned <- metadata_new[match(gsva_ids, metadata_new$SampleID), ]
  stopifnot(all(gsva_ids == metadata_aligned$SampleID))
  
  gsva_df <- data.frame(SampleID = gsva_ids,
                        GSVA_score = gsva_vec,
                        Phenotype = metadata_aligned[[phenotype_col]])
  
  # applying threshold for classification.
  pred <- ifelse(gsva_df$GSVA_score > threshold, label_one, label_two)
  
  # confusion matrices.
  cm <- confusionMatrix(factor(pred, levels = c(label_one, label_two)),
                        factor(gsva_df$Phenotype, levels = c(label_one, label_two)))
  
  # 
  roc_ext <- roc(response = gsva_df$Phenotype,
                 predictor = gsva_df$GSVA_score,
                 levels = c(label_two, label_one),
                 direction = "<")
  
  # result.
  results <- list(
    auc = as.numeric(auc(roc_ext)),
    threshold_used = threshold,
    metrics = data.frame(
      Accuracy  = cm$overall["Accuracy"],
      Precision = cm$byClass["Precision"],
      Recall    = cm$byClass["Recall"],
      F1        = cm$byClass["F1"]
    ),
    confusion = cm$table,
    gsva_scores = gsva_df
  )
  
  return(results)
}

res2 <- evaluate_gsva_external(expr_new = ackerman_NB_log2, metadata_new = ackerman_metadata,
                               geneset_tmmUp = best_geneset1, geneset_tmmDown = best_geneset2,
                               phenotype_col = "TMM_Case",
                               label_one = "NO_TMM", label_two = "TMM",
                               threshold = 0.6316002,
                               kcdf = "Gaussian")

