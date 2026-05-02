#!/usr/bin/env Rscript
# Model Performance Evaluation with Bootstrap CI

library(data.table)
library(pROC)
library(dplyr)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript evaluate_model.R <predictions_dir> <output_file>")
}

pred_dir <- args[1]
output_file <- args[2]

# Find optimal cutoff
find_optimal_cutoff <- function(y_true, y_pred) {
  roc_obj <- roc(y_true, y_pred, quiet = TRUE)
  coords_obj <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
  return(coords_obj$threshold)
}

# Calculate metrics
calculate_metrics <- function(y_true, y_pred, cutoff) {
  y_pred_binary <- ifelse(y_pred >= cutoff, 1, 0)
  cm <- table(Predicted = y_pred_binary, Actual = y_true)
  
  if (nrow(cm) == 2 && ncol(cm) == 2) {
    tn <- cm[1,1]; fp <- cm[1,2]; fn <- cm[2,1]; tp <- cm[2,2]
  } else {
    tn <- fp <- fn <- tp <- 0
  }
  
  acc <- (tp + tn) / (tp + tn + fp + fn)
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), 0)
  auc <- as.numeric(roc(y_true, y_pred, quiet = TRUE)$auc)
  
  return(c(auc, acc, sens, spec))
}

# Bootstrap evaluation
bootstrap_eval <- function(y_true, y_pred, n_iter = 1000) {
  cutoff <- find_optimal_cutoff(y_true, y_pred)
  n <- length(y_true)
  results <- matrix(0, nrow = n_iter, ncol = 4)
  
  set.seed(123)
  for (i in 1:n_iter) {
    boot_idx <- sample(1:n, n, replace = TRUE)
    results[i, ] <- calculate_metrics(y_true[boot_idx], y_pred[boot_idx], cutoff)
  }
  
  median_vals <- apply(results, 2, median)
  ci_lower <- apply(results, 2, quantile, 0.025)
  ci_upper <- apply(results, 2, quantile, 0.975)
  
  formatted <- paste0(
    sprintf("%.3f", median_vals), " [",
    sprintf("%.3f", ci_lower), " - ",
    sprintf("%.3f", ci_upper), "]"
  )
  
  names(formatted) <- c('AUC', 'Accuracy', 'Sensitivity', 'Specificity')
  return(formatted)
}

# Process all prediction files
pred_files <- list.files(pred_dir, pattern = "^predictions_.*\\.csv$", full.names = TRUE)
cat("Found", length(pred_files), "prediction files\n")

results <- list()

for (pred_file in pred_files) {
  disease <- gsub("predictions_(.*)\\.csv", "\\1", basename(pred_file))
  cat("Processing:", disease, "\n")
  
  tryCatch({
    pred_data <- fread(pred_file)
    
    if (!"true_label" %in% colnames(pred_data)) {
      cat("  Warning: No true_label column, skipping\n")
      next
    }
    
    # Find prediction columns (exclude ID and label)
    pred_cols <- setdiff(colnames(pred_data), c("sample_id", "true_label", "ID"))
    
    for (pred_col in pred_cols) {
      valid_mask <- complete.cases(pred_data[[pred_col]], pred_data$true_label)
      
      if (sum(valid_mask) >= 50) {
        y_true <- pred_data$true_label[valid_mask]
        y_pred <- pred_data[[pred_col]][valid_mask]
        
        metrics <- bootstrap_eval(y_true, y_pred)
        
        results[[length(results) + 1]] <- data.frame(
          disease = disease,
          model = pred_col,
          AUC = metrics["AUC"],
          Accuracy = metrics["Accuracy"],
          Sensitivity = metrics["Sensitivity"],
          Specificity = metrics["Specificity"],
          stringsAsFactors = FALSE
        )
        
        cat("  ✓", pred_col, "- AUC:", strsplit(metrics["AUC"], " ")[[1]][1], "\n")
      }
    }
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
  })
}

# Save results 
if (length(results) > 0) {
  final_df <- do.call(rbind, results)
  fwrite(final_df, output_file)
  cat("\n✓ Results saved to:", output_file, "\n")
  cat("Total records:", nrow(final_df), "\n")
} else {
  cat("\nWarning: No results generated\n")
}

# DeLong检验
test_result <- roc.test(roc_1, roc_2, method = "delong")