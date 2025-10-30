
## finding the best gene signature for NO_TMM using AUC -- greedy forward selection.

library(pROC)
library(caret)
library(GSVA)

run_exhaustive_forward_auc <- function(expr_matrix, metadata, candidate_genes, 
                                       phenotype_col,
                                       label_one, label_two, # labels would be NO_TMM or TMM in our case.
                                       max_genes,
                                       pivot_gene)  {  # gene giving the best results from linear regression test.
  expr_matrix <- as.matrix(expr_matrix)
  stopifnot("SampleID" %in% colnames(metadata))
  
  # aligning samples.
  group_labels <- metadata[[phenotype_col]]
  one_samples  <- metadata$SampleID[group_labels == label_one]
  two_samples  <- metadata$SampleID[group_labels == label_two]
  keep_samples <- sort(intersect(colnames(expr_matrix), c(one_samples, two_samples)))
  
  binary_labels <- setNames(ifelse(keep_samples %in% one_samples, 1, 0), keep_samples) # for auc.
  
  
  # helper function: GSVA score for a gene set on aligned samples
  .gsva_scores <- function(genes) {
    
    
    param <- gsvaParam(
      exprData = expr_matrix[, keep_samples, drop = FALSE],
      geneSets = list(sig = genes),
      kcdf     = "Gaussian" # for our log counts data.
    )
    mat <- gsva(param, verbose = FALSE)
    as.numeric(mat["sig", ])
  }
  
  selected_genes <- intersect(pivot_gene, rownames(expr_matrix))
  
  seed_scores <- .gsva_scores(selected_genes)
  if (is.null(seed_scores)) {
    # falling back: pick the best single-gene seed from candidates with non-zero variance.
    pool <- intersect(candidate_genes, rownames(expr_matrix))
    # drop zero-variance candidates
    zv <- apply(expr_matrix[pool, keep_samples, drop=FALSE], 1, function(x) sd(x, na.rm=TRUE) == 0)
    pool <- pool[!zv]
    if (length(pool) == 0) 
      stop("No candidate genes with non-zero variance in kept samples.")
    
    # picking the one with highest AUC.
    best_auc <- -Inf
    best_gene <- NULL
    best_scores <- NULL
    
    for (g in pool) {
      sc <- .gsva_scores(g)
      if (is.null(sc)) next
      names(sc) <- keep_samples
      a <- as.numeric(auc(roc(binary_labels, sc, quiet = TRUE)))
      if (!is.na(a) && a > best_auc) {
        best_auc <- a; 
        best_gene <- g; 
        best_scores <- sc
      }
    }
    if (is.null(best_gene)) 
      stop("Unable to initialize a seed gene for GSVA.")
    
    selected_genes <- best_gene
    scores0 <- best_scores
    auc0 <- best_auc
    message("Best combination so far: ", best_gene, " | AUC = ", round(auc0, 4))
  } 
  
  else {
    names(seed_scores) <- keep_samples
    auc0 <- as.numeric(auc(roc(binary_labels, seed_scores, quiet = TRUE)))
    scores0 <- seed_scores
    message("Baseline (", pivot_gene, "): AUC = ", round(auc0, 4))
  }
  
  best_auc <- auc0
  results <- list()
  results[[1]] <- list(genes = selected_genes, auc = best_auc, scores = scores0)
  
  # remaining genes.
  remaining_genes <- setdiff(intersect(candidate_genes, rownames(expr_matrix)), selected_genes)
  
  # forward selection.
  while (length(selected_genes) < max_genes && length(remaining_genes) > 0) {
    best_iteration_auc <- -Inf
    best_gene <- NULL
    best_result <- NULL
    
    for (gene in remaining_genes) {
      test_genes <- c(selected_genes, gene)
      sc <- .gsva_scores(test_genes)
      if (is.null(sc)) 
        next
      names(sc) <- keep_samples
      
      roc_obj <- roc(binary_labels, sc, quiet = TRUE)
      current_auc <- as.numeric(auc(roc_obj))
      
      if (!is.na(current_auc) && current_auc > best_iteration_auc) {
        best_iteration_auc <- current_auc
        best_gene <- gene
        best_result <- list(
          genes  = test_genes,
          auc    = current_auc,
          roc    = roc_obj,
          scores = sc
        )
      }
    }
    
    if (is.null(best_gene) || best_iteration_auc <= best_auc) {
      message("Stopping: No AUC improvement at size ", length(selected_genes) + 1)
      break
    }
    
    selected_genes <- c(selected_genes, best_gene)
    best_auc <- best_iteration_auc
    results[[length(selected_genes)]] <- best_result
    remaining_genes <- setdiff(remaining_genes, best_gene)
    
    message("Step ", length(selected_genes), ": Added ", best_gene, " | AUC = ", round(best_auc, 4))
  }
  
  list(final_genes = selected_genes, final_auc = best_auc, results = results)
}



#########################################################################################################

## function for running k-fold validation.

run_kfold_auc <- function(expr_matrix, metadata, candidate_genes,
                          phenotype_col, label_one, label_two,
                          max_genes, k,
                          pivot_gene) {
  set.seed(123)
  
  # aligning samples across expression and metadata.
  expr_matrix <- as.matrix(expr_matrix)
  stopifnot("SampleID" %in% colnames(metadata))
  
  common_samples <- intersect(colnames(expr_matrix), metadata$SampleID)
  if (length(common_samples) < 4) stop("Not enough overlapping samples between expr_matrix and metadata.")
  
  expr_matrix <- expr_matrix[, common_samples, drop = FALSE]
  metadata    <- metadata[match(common_samples, metadata$SampleID), ]
  
  # Building y (factor with two levels).
  y <- metadata[[phenotype_col]]
  if (anyNA(y)) {
    keep <- !is.na(y)
    expr_matrix <- expr_matrix[, keep, drop = FALSE]
    metadata    <- metadata[keep, , drop = FALSE]
    y <- metadata[[phenotype_col]]
  }
  y <- factor(y)
  if (nlevels(y) < 2) 
    stop("phenotype_col has < 2 classes after alignment.")
  if (!all(c(label_one, label_two) %in% levels(y))) {
    stop("label_one/label_two not found in phenotype_col levels after alignment.")
  }
  
  # Adjust k if it's larger than the minority class count
  class_counts <- table(y)
  max_k <- min(class_counts)
  if (k > max_k) {
    message("Requested k (", k, ") > smallest class size (", max_k, "). Using k = ", max_k, ".")
    k <- as.integer(max_k)
  }
  if (k < 2) stop("k must be >= 2 and <= min class size.")
  
  # Create stratified folds on aligned y
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  if (length(folds) == 0) stop("createFolds returned 0 folds. Check phenotype_col and k.")
  
  fold_results <- vector("list", length(folds))
  fold_aucs    <- rep(NA_real_, length(folds))
  
  for (i in seq_along(folds)) {
    message("Fold: ", i)
    
    test_idx     <- folds[[i]]
    test_samples <- metadata$SampleID[test_idx]
    train_samples <- setdiff(metadata$SampleID, test_samples)
    
    # Train split
    expr_train     <- expr_matrix[, train_samples, drop = FALSE]
    metadata_train <- metadata[match(train_samples, metadata$SampleID), , drop = FALSE]
    
    # Running forward GSVA selection on the TRAIN ONLY.
    res <- run_exhaustive_forward_auc(
      expr_matrix   = expr_train,
      metadata      = metadata_train,
      candidate_genes = candidate_genes,
      phenotype_col = phenotype_col,
      label_one     = label_one,
      label_two     = label_two,
      max_genes     = max_genes,
      pivot_gene     = pivot_gene
    )
    selected_genes <- res$final_genes
    
    # Test split
    expr_test     <- expr_matrix[, test_samples, drop = FALSE]
    metadata_test <- metadata[match(test_samples, metadata$SampleID), , drop = FALSE]
    
    # Computing GSVA scores on TEST.
    param <- gsvaParam(
      exprData = expr_test,
      geneSets = list(sig = selected_genes),
      kcdf     = "Gaussian",
      minSize  = 1L
    )
    gsva_mat <- gsva(param, verbose = FALSE)
    scores   <- as.numeric(gsva_mat["sig", ])
    names(scores) <- colnames(expr_test)
    
    # Labels on TEST
    y_test <- ifelse(metadata_test[[phenotype_col]] == label_one, 1, 0)
    names(y_test) <- metadata_test$SampleID
    
    common <- intersect(names(scores), names(y_test))
    if (length(common) >= 2 && length(unique(y_test[common])) == 2) {
      roc_obj   <- roc(y_test[common], scores[common], quiet = TRUE)
      fold_auc  <- as.numeric(auc(roc_obj))
      fold_aucs[i] <- fold_auc
    } else {
      fold_aucs[i] <- NA_real_
      message("  (Insufficient class variation in test fold; AUC set to NA)")
    }
    
    fold_results[[i]] <- list(
      selected_genes = selected_genes,
      fold_auc       = fold_aucs[i]
    )
  }
  
  list(
    fold_results = fold_results,
    aucs         = fold_aucs,
    mean_auc     = mean(fold_aucs, na.rm = TRUE),
    sd_auc       = sd(fold_aucs,  na.rm = TRUE)
  )
}


####################









