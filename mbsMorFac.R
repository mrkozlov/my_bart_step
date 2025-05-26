# File: my_bart_step_updated_with_summary.R
# Description: Updated BART modeling script with Final Summary section for small sample sizes with overfitting analysis, variable importance, and spatial analysis.

library(embarcadero)
library(dbarts)
library(pROC)
library(ggplot2)
library(spdep)
library(gridExtra)
library(reshape2)

options(error = function() {
  if (sink.number() > 0) sink()
  while (dev.cur() > 1) dev.off()
  closeAllConnections()
})

calculate_metrics <- function(true_labels, predictions, threshold) {
  pred_labels <- ifelse(predictions >= threshold, 1, 0)
  cm <- table(True = true_labels, Predicted = pred_labels)
  tn <- ifelse("0" %in% rownames(cm) && "0" %in% colnames(cm), cm["0", "0"], 0)
  fp <- ifelse("0" %in% rownames(cm) && "1" %in% colnames(cm), cm["0", "1"], 0)
  fn <- ifelse("1" %in% rownames(cm) && "0" %in% colnames(cm), cm["1", "0"], 0)
  tp <- ifelse("1" %in% rownames(cm) && "1" %in% colnames(cm), cm["1", "1"], 0)
  
  fnr <- ifelse((fn + tp) > 0, fn / (fn + tp), 0)
  fpr <- ifelse((fp + tn) > 0, fp / (fp + tn), 0)
  tss <- (tp / (tp + fn)) + (tn / (tn + fp)) - 1
  
  return(list(tss = tss, fnr = fnr, fpr = fpr, cm = cm))
}

find_threshold_for_fnr <- function(predictions, true_labels, target_fnr) {
  n_pos <- sum(true_labels == 1, na.rm = TRUE)
  if (n_pos == 0) {
    warning("No positive cases in the data. Returning NA for threshold.")
    return(list(threshold = NA, metrics = NA))
  }
  
  max_fnr <- 1 / n_pos
  if (target_fnr > 0 && target_fnr < 1) {
    target_fnr <- max(0.04, min(0.2, max_fnr))
    cat("Adjusted target FNR for sample size n =", n_pos, "to", target_fnr, "\n")
  }
  
  if (target_fnr == 0) {
    threshold <- min(predictions[true_labels == 1], na.rm = TRUE)
    metrics <- calculate_metrics(true_labels, predictions, threshold)
    return(list(threshold = threshold, metrics = metrics))
  }
  
  all_preds <- sort(unique(predictions), decreasing = FALSE)
  prev_fnr <- 0
  prev_threshold <- NA
  
  for (i in seq_along(all_preds)) {
    threshold <- all_preds[i]
    metrics <- calculate_metrics(true_labels, predictions, threshold)
    
    if (!is.na(metrics$fnr) && metrics$fnr >= target_fnr) {
      if (i > 1 && prev_fnr == 0 && metrics$fnr > 0) {
        threshold <- (prev_threshold + threshold) / 2
        metrics <- calculate_metrics(true_labels, predictions, threshold)
      }
      return(list(threshold = threshold, metrics = metrics))
    }
    prev_fnr <- metrics$fnr
    prev_threshold <- threshold
  }
  
  fnr0_threshold <- min(predictions[true_labels == 1], na.rm = TRUE)
  next_threshold_idx <- which(all_preds > fnr0_threshold)[1]
  if (!is.na(next_threshold_idx)) {
    threshold <- all_preds[next_threshold_idx]
    metrics <- calculate_metrics(true_labels, predictions, threshold)
    cat("Forcing threshold above FNR=0 to achieve FNR > 0: ", threshold, "\n")
    return(list(threshold = threshold, metrics = metrics))
  } else {
    threshold <- fnr0_threshold + 0.01
    metrics <- calculate_metrics(true_labels, predictions, threshold)
    cat("No higher threshold available, incrementing by 0.01 to: ", threshold, "\n")
    return(list(threshold = threshold, metrics = metrics))
  }
}

find_optimal_k <- function(residuals, coords, max_k = NULL, min_points = 3) {
  n <- length(residuals)
  if (is.null(max_k)) {
    max_k <- max(min_points, min(floor(n / 3), 10))
  }
  
  k_range <- 1:max_k
  best_k <- NULL
  best_p_value <- 1
  moran_results <- list()
  
  cat("\n=== Moran's I for different k values ===\n")
  for (k in k_range) {
    nb <- tryCatch({
      knn2nb(knearneigh(coords, k = k))
    }, error = function(e) {
      warning("Error in knearneigh for k = ", k, ": ", e$message)
      return(NULL)
    })
    if (is.null(nb)) next
    
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(residuals, listw, zero.policy = TRUE)
    p_value <- moran_test$p.value
    
    cat("k =", k, ":\n")
    cat("Moran I =", moran_test$estimate[1], ", p-value =", p_value, "\n\n")
    
    moran_results[[as.character(k)]] <- moran_test
    if (p_value < best_p_value) {
      best_p_value <- p_value
      best_k <- k
    }
  }
  
  if (is.null(best_k)) {
    warning("No valid k found. Using default k = 1.")
    best_k <- 1
  }
  
  return(list(optimal_k = best_k, p_value = best_p_value, results = moran_results))
}

cross_moran <- function(residuals_train, residuals_test, coords_train, coords_test, k) {
  nb_train_to_test <- knn2nb(knearneigh(coords_train, k = k, coords_test))
  listw_train_to_test <- nb2listw(nb_train_to_test, style = "W", zero.policy = TRUE)
  cross_moran_test <- moran.test(residuals_train, listw_train_to_test, zero.policy = TRUE)
  return(cross_moran_test)
}

my_bart_step <- function(x.data, y.data, test_data, 
                         selected_vars = setdiff(colnames(x.data), "nf"),
                         tree.step = 5, iter.plot = 3, iter.step = 2,
                         full = TRUE, quiet = FALSE,
                         num_trees = 200, num_burn_in = 250,
                         num_iterations_after_burn_in = 1000,
                         target_fnr = 0.04,
                         glm_family = "binomial",
                         auc_weight = 0.5, ks_weight = 0.5, ...) {
  
  sink("bart_report.txt")
  on.exit(if (sink.number() > 0) sink(), add = TRUE)
  
  # Input Data Checks
  cat("Checking for NA in input data:\n")
  cat("NA in y.data:", sum(is.na(y.data)), "\n")
  cat("NA in test_data$nf:", sum(is.na(test_data$nf)), "\n")
  cat("NA in x.data:", sum(is.na(x.data[, selected_vars])), "\n")
  
  if (sum(is.na(y.data)) > 0 || sum(is.na(test_data$nf)) > 0 || sum(is.na(x.data[, selected_vars])) > 0) {
    warning("Input data contains NA values. Consider handling them before proceeding.")
  }
  
  cat("Checking training data:\n")
  cat("Number of observations:", length(y.data), "\n")
  print(table(y.data))
  cat("Positive (1):", sum(y.data == 1), ", Negative (0):", sum(y.data == 0), ", NA:", sum(is.na(y.data)), "\n")
  
  cat("\nChecking test data:\n")
  cat("Number of observations:", nrow(test_data), "\n")
  print(table(test_data$nf))
  cat("Positive (1):", sum(test_data$nf == 1), ", Negative (0):", sum(test_data$nf == 0), ", NA:", sum(is.na(test_data$nf)), "\n")
  
  cat("\nChecking selected_vars:\n")
  cat("Variables in x.data:", paste(colnames(x.data), collapse = " "), "\n")
  cat("Selected variables:", paste(selected_vars, collapse = " "), "\n")
  if (!all(selected_vars %in% colnames(x.data))) {
    stop("Some variables from selected_vars are missing in x.data")
  }
  
  # Define do_spatial_analysis
  do_spatial_analysis <- all(c("x", "y") %in% colnames(x.data)) && all(c("x", "y") %in% colnames(test_data))
  
  # Run bart.step
  cat("\n=== Running bart.step ===\n")
  pdf("bart_step_diagnostic_plots.pdf", width = 8, height = 6)
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  output <- capture.output({
    bart_model <- suppressWarnings(bart.step(
      x.data = x.data[, selected_vars, drop = FALSE], 
      y.data = y.data,
      tree.step = tree.step, 
      iter.plot = iter.plot, 
      iter.step = iter.step,
      full = full, 
      quiet = quiet, 
      ...
    ))
  })
  
  cat("Diagnostic plots saved to bart_step_diagnostic_plots.pdf\n")
  cat("Captured output from bart.step:\n")
  cat(paste(output, collapse = "\n"), "\n")
  
  # Extract Final Variables
  cat("\n=== Debugging bart.step output ===\n")
  final_vars <- NULL
  if (any(grepl("Final recommended variable list", output))) {
    cat("Found 'Final recommended variable list' in output.\n")
    final_vars_idx <- grep("Final recommended variable list", output)
    if ((final_vars_idx + 1) <= length(output)) {
      line <- output[final_vars_idx + 1]
      cat("Raw line after 'Final recommended variable list': '", line, "'\n")
      line <- gsub("\\[1\\]\\s*", "", line)
      line <- trimws(line)
      cat("Processed line: '", line, "'\n")
      final_vars <- unlist(strsplit(line, "\\s+"))
      final_vars <- final_vars[final_vars != "" & final_vars != "nf"]
      cat("Extracted final_vars: ", paste(final_vars, collapse = ", "), "\n")
    } else {
      cat("No line found after 'Final recommended variable list'.\n")
    }
  } else {
    cat("Could not find 'Final recommended variable list' in output. Possible reasons:\n")
    cat("- bart.step did not complete variable selection.\n")
    cat("- Output format changed in this version of embarcadero.\n")
    cat("- Insufficient iterations or data issues.\n")
  }
  
  if (is.null(final_vars) || length(final_vars) == 0) {
    warning("Could not find variables after 'Final recommended variable list' in output. Using selected_vars.")
    final_vars <- selected_vars
  }
  cat("Variables used (bart.step):", paste(final_vars, collapse = ", "), "\n")
  
  # Extract Thresholds and Metrics from Output
  bart_thr_tss_idx <- grep("Recommended threshold", output)
  bart_fnr_idx <- grep("Resulting type I error rate", output)
  bart_fpr_idx <- grep("Resulting type II error rate", output)
  
  bart_thr_tss <- if (length(bart_thr_tss_idx) > 0) as.numeric(gsub("Cutoff =\\s*", "", output[bart_thr_tss_idx + 1])) else NA
  bart_fnr_from_output <- if (length(bart_fnr_idx) > 0) as.numeric(gsub(".*:\\s*", "", output[bart_fnr_idx + 1])) else NA
  bart_fpr_from_output <- if (length(bart_fpr_idx) > 0) as.numeric(gsub(".*:\\s*", "", output[bart_fpr_idx + 1])) else NA
  
  bart_train_pred <- fitted(bart_model)
  bart_train_pred <- pmin(pmax(bart_train_pred, 0), 1)
  
  # Calculate Metrics for bart.step (Train)
  roc_bart_train <- roc(y.data, bart_train_pred, quiet = TRUE)
  auc_bart_train <- auc(roc_bart_train)
  bart_fnr0_result <- find_threshold_for_fnr(bart_train_pred, y.data, 0)
  bart_thr_fnr0 <- bart_fnr0_result$threshold
  bart_metrics_fnr0 <- bart_fnr0_result$metrics
  bart_fnr05_result <- find_threshold_for_fnr(bart_train_pred, y.data, target_fnr)
  bart_thr_fnr05 <- bart_fnr05_result$threshold
  bart_metrics_fnr05 <- bart_fnr05_result$metrics
  bart_coords <- coords(roc_bart_train, "best", ret = "threshold")
  if (is.na(bart_thr_tss)) {
    bart_thr_tss <- bart_coords$threshold[1]
  }
  bart_metrics_tss <- calculate_metrics(y.data, bart_train_pred, bart_thr_tss)
  
  if (is.na(bart_fpr_from_output) || is.na(bart_fnr_from_output)) {
    bart_metrics_tss_manual <- calculate_metrics(y.data, bart_train_pred, bart_thr_tss)
    bart_fpr_from_output <- bart_metrics_tss_manual$fpr
    bart_fnr_from_output <- bart_metrics_tss_manual$fnr
    bart_metrics_tss <- list(tss = 1 - bart_fnr_from_output + (1 - bart_fpr_from_output) - 1,
                             fnr = bart_fnr_from_output, fpr = bart_fpr_from_output)
  }
  
  train_predictors <- x.data[, final_vars, drop = FALSE]
  test_predictors <- test_data[, final_vars, drop = FALSE]
  cat("\nTrain predictors columns:", paste(colnames(train_predictors), collapse = ", "), "\n")
  cat("Test predictors columns:", paste(colnames(test_predictors), collapse = ", "), "\n")
  
  # Run wbart
  cat("\n=== Running BART::wbart ===\n")
  wbart_model <- bart(x.train = train_predictors, y.train = y.data, x.test = test_predictors,
                      ntree = num_trees, ndpost = num_iterations_after_burn_in, nskip = num_burn_in,
                      keeptrees = TRUE, verbose = FALSE)
  wbart_train_pred <- plogis(colMeans(wbart_model$yhat.train))
  wbart_train_pred <- pmin(pmax(wbart_train_pred, 0), 1)
  wbart_test_pred <- plogis(colMeans(wbart_model$yhat.test))
  wbart_test_pred <- pmin(pmax(wbart_test_pred, 0), 1)
  cat("Unique wbart_test_pred:", paste(sort(unique(wbart_test_pred)), collapse = " "), "\n")
  
  # Calculate Metrics for wbart (Train)
  roc_wbart_train <- roc(y.data, wbart_train_pred, quiet = TRUE)
  auc_wbart_train <- auc(roc_wbart_train)
  wbart_fnr0_train_result <- find_threshold_for_fnr(wbart_train_pred, y.data, 0)
  wbart_thr_fnr0_train <- wbart_fnr0_train_result$threshold
  wbart_metrics_fnr0_train <- wbart_fnr0_train_result$metrics
  wbart_fnr05_train_result <- find_threshold_for_fnr(wbart_train_pred, y.data, target_fnr)
  wbart_thr_fnr05_train <- wbart_fnr05_train_result$threshold
  wbart_metrics_fnr05_train <- wbart_fnr05_train_result$metrics
  wbart_coords_train <- coords(roc_wbart_train, "best", ret = "threshold")
  wbart_thr_tss_train <- wbart_coords_train$threshold[1]
  wbart_metrics_tss_train <- calculate_metrics(y.data, wbart_train_pred, wbart_thr_tss_train)
  
  # Calculate Metrics for wbart (Test)
  roc_wbart_test <- roc(test_data$nf, wbart_test_pred, quiet = TRUE)
  auc_wbart_test <- auc(roc_wbart_test)
  wbart_fnr0_result <- find_threshold_for_fnr(wbart_test_pred, test_data$nf, 0)
  wbart_thr_fnr0 <- wbart_fnr0_result$threshold
  wbart_metrics_fnr0 <- wbart_fnr0_result$metrics
  wbart_fnr05_result <- find_threshold_for_fnr(wbart_test_pred, test_data$nf, target_fnr)
  wbart_thr_fnr05 <- wbart_fnr05_result$threshold
  wbart_metrics_fnr05 <- wbart_fnr05_result$metrics
  wbart_coords_test <- coords(roc_wbart_test, "best", ret = "threshold")
  wbart_thr_tss <- wbart_coords_test$threshold[1]
  wbart_metrics_tss <- calculate_metrics(test_data$nf, wbart_test_pred, wbart_thr_tss)
  
  # Run dbarts
  cat("\n=== Running dbarts::bart ===\n")
  dbart_model <- dbarts::bart(x.train = train_predictors, y.train = y.data, x.test = test_predictors,
                              ntree = num_trees, ndpost = num_iterations_after_burn_in, nskip = num_burn_in,
                              keeptrees = TRUE, verbose = FALSE)
  dbart_train_pred <- plogis(colMeans(dbart_model$yhat.train))
  dbart_train_pred <- pmin(pmax(dbart_train_pred, 0), 1)
  dbart_test_pred <- plogis(colMeans(dbart_model$yhat.test))
  dbart_test_pred <- pmin(pmax(dbart_test_pred, 0), 1)
  cat("Unique dbart_test_pred:", paste(sort(unique(dbart_test_pred)), collapse = " "), "\n")
  
  # Calculate Metrics for dbart (Train)
  roc_dbart_train <- roc(y.data, dbart_train_pred, quiet = TRUE)
  auc_dbart_train <- auc(roc_dbart_train)
  dbart_fnr0_train_result <- find_threshold_for_fnr(dbart_train_pred, y.data, 0)
  dbart_thr_fnr0_train <- dbart_fnr0_train_result$threshold
  dbart_metrics_fnr0_train <- dbart_fnr0_train_result$metrics
  dbart_fnr05_train_result <- find_threshold_for_fnr(dbart_train_pred, y.data, target_fnr)
  dbart_thr_fnr05_train <- dbart_fnr05_train_result$threshold
  dbart_metrics_fnr05_train <- dbart_fnr05_train_result$metrics
  dbart_coords_train <- coords(roc_dbart_train, "best", ret = "threshold")
  dbart_thr_tss_train <- dbart_coords_train$threshold[1]
  dbart_metrics_tss_train <- calculate_metrics(y.data, dbart_train_pred, dbart_thr_tss_train)
  
  # Calculate Metrics for dbart (Test)
  roc_dbart_test <- roc(test_data$nf, dbart_test_pred, quiet = TRUE)
  auc_dbart_test <- auc(roc_dbart_test)
  dbart_fnr0_result <- find_threshold_for_fnr(dbart_test_pred, test_data$nf, 0)
  dbart_thr_fnr0 <- dbart_fnr0_result$threshold
  dbart_metrics_fnr0 <- dbart_fnr0_result$metrics
  dbart_fnr05_result <- find_threshold_for_fnr(dbart_test_pred, test_data$nf, target_fnr)
  dbart_thr_fnr05 <- dbart_fnr05_result$threshold
  dbart_metrics_fnr05 <- dbart_fnr05_result$metrics
  dbart_coords_test <- coords(roc_dbart_test, "best", ret = "threshold")
  dbart_thr_tss <- dbart_coords_test$threshold[1]
  dbart_metrics_tss <- calculate_metrics(test_data$nf, dbart_test_pred, dbart_thr_tss)
  
  # Calibration for bart.step (Test)
  if (glm_family == "binomial") {
    glm_calib <- tryCatch({
      glm(y.data ~ wbart_train_pred + dbart_train_pred, family = "binomial")
    }, warning = function(w) {
      message("Warning in GLM with family = 'binomial': ", w$message)
      message("Switching to quasibinomial family.")
      glm(y.data ~ wbart_train_pred + dbart_train_pred, family = "quasibinomial")
    })
  } else {
    glm_calib <- glm(y.data ~ wbart_train_pred + dbart_train_pred, family = "quasibinomial")
  }
  bart_test_pred_glm <- predict(glm_calib, newdata = data.frame(wbart_train_pred = wbart_test_pred, dbart_train_pred = dbart_test_pred), type = "response")
  bart_test_pred_glm <- pmin(pmax(bart_test_pred_glm, 0), 1)
  cat("Unique bart_test_pred_glm:", paste(sort(unique(bart_test_pred_glm)), collapse = " "), "\n")
  
  lm_calib <- lm(bart_train_pred ~ wbart_train_pred + dbart_train_pred)
  bart_test_pred_lm <- predict(lm_calib, newdata = data.frame(wbart_train_pred = wbart_test_pred, dbart_train_pred = dbart_test_pred))
  bart_test_pred_lm <- pmin(pmax(bart_test_pred_lm, 0), 1)
  cat("Unique bart_test_pred_lm:", paste(sort(unique(bart_test_pred_lm)), collapse = " "), "\n")
  
  bart_test_pred_quant <- dbart_test_pred
  
  # Calculate ROC and AUC for Calibrated Models
  roc_bart_test_quant <- roc(test_data$nf, bart_test_pred_quant, quiet = TRUE)
  auc_bart_test_quant <- auc(roc_bart_test_quant)
  bart_fnr0_test_quant_result <- find_threshold_for_fnr(bart_test_pred_quant, test_data$nf, 0)
  bart_thr_fnr0_test_quant <- bart_fnr0_test_quant_result$threshold
  bart_metrics_fnr0_test_quant <- bart_fnr0_test_quant_result$metrics
  bart_fnr05_test_quant_result <- find_threshold_for_fnr(bart_test_pred_quant, test_data$nf, target_fnr)
  bart_thr_fnr05_test_quant <- bart_fnr05_test_quant_result$threshold
  bart_metrics_fnr05_test_quant <- bart_fnr05_test_quant_result$metrics
  bart_coords_test_quant <- coords(roc_bart_test_quant, "best", ret = "threshold")
  bart_thr_tss_test_quant <- bart_coords_test_quant$threshold[1]
  bart_metrics_tss_test_quant <- calculate_metrics(test_data$nf, bart_test_pred_quant, bart_thr_tss_test_quant)
  
  roc_bart_test_glm <- roc(test_data$nf, bart_test_pred_glm, quiet = TRUE)
  auc_bart_test_glm <- auc(roc_bart_test_glm)
  bart_fnr0_test_glm_result <- find_threshold_for_fnr(bart_test_pred_glm, test_data$nf, 0)
  bart_thr_fnr0_test_glm <- bart_fnr0_test_glm_result$threshold
  bart_metrics_fnr0_test_glm <- bart_fnr0_test_glm_result$metrics
  bart_fnr05_test_glm_result <- find_threshold_for_fnr(bart_test_pred_glm, test_data$nf, target_fnr)
  bart_thr_fnr05_test_glm <- bart_fnr05_test_glm_result$threshold
  bart_metrics_fnr05_test_glm <- bart_fnr05_test_glm_result$metrics
  bart_coords_test_glm <- coords(roc_bart_test_glm, "best", ret = "threshold")
  bart_thr_tss_test_glm <- bart_coords_test_glm$threshold[1]
  bart_metrics_tss_test_glm <- calculate_metrics(test_data$nf, bart_test_pred_glm, bart_thr_tss_test_glm)
  
  roc_bart_test_lm <- roc(test_data$nf, bart_test_pred_lm, quiet = TRUE)
  auc_bart_test_lm <- auc(roc_bart_test_lm)
  bart_fnr0_test_lm_result <- find_threshold_for_fnr(bart_test_pred_lm, test_data$nf, 0)
  bart_thr_fnr0_test_lm <- bart_fnr0_test_lm_result$threshold
  bart_metrics_fnr0_test_lm <- bart_fnr0_test_lm_result$metrics
  bart_fnr05_test_lm_result <- find_threshold_for_fnr(bart_test_pred_lm, test_data$nf, target_fnr)
  bart_thr_fnr05_test_lm <- bart_fnr05_test_lm_result$threshold
  bart_metrics_fnr05_test_lm <- bart_fnr05_test_lm_result$metrics
  bart_coords_test_lm <- coords(roc_bart_test_lm, "best", ret = "threshold")
  bart_thr_tss_test_lm <- bart_coords_test_lm$threshold[1]
  bart_metrics_tss_test_lm <- calculate_metrics(test_data$nf, bart_test_pred_lm, bart_thr_tss_test_lm)
  
  # Check Predictions for NA
  cat("\nChecking predictions for NA:\n")
  cat("NA in bart_train_pred:", sum(is.na(bart_train_pred)), "\n")
  cat("NA in wbart_train_pred:", sum(is.na(wbart_train_pred)), "\n")
  cat("NA in dbart_train_pred:", sum(is.na(dbart_train_pred)), "\n")
  cat("NA in wbart_test_pred:", sum(is.na(wbart_test_pred)), "\n")
  cat("NA in dbart_test_pred:", sum(is.na(dbart_test_pred)), "\n")
  cat("NA in bart_test_pred_quant:", sum(is.na(bart_test_pred_quant)), "\n")
  cat("NA in bart_test_pred_glm:", sum(is.na(bart_test_pred_glm)), "\n")
  cat("NA in bart_test_pred_lm:", sum(is.na(bart_test_pred_lm)), "\n")
  
  # Output Metrics for All Models
  cat("\n=== Model Metrics ===\n")
  
  cat("\n=== bart.step (Train) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_bart_train, bart_thr_fnr0, bart_metrics_fnr0$tss, bart_metrics_fnr0$fnr, bart_metrics_fnr0$fpr,
              target_fnr, bart_thr_fnr05, bart_metrics_fnr05$tss, bart_metrics_fnr05$fnr, bart_metrics_fnr05$fpr,
              bart_thr_tss, bart_metrics_tss$tss, bart_metrics_tss$fnr, bart_metrics_tss$fpr))
  
  cat("\n=== wbart (Train) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_wbart_train, wbart_thr_fnr0_train, wbart_metrics_fnr0_train$tss, wbart_metrics_fnr0_train$fnr, wbart_metrics_fnr0_train$fpr,
              target_fnr, wbart_thr_fnr05_train, wbart_metrics_fnr05_train$tss, wbart_metrics_fnr05_train$fnr, wbart_metrics_fnr05_train$fpr,
              wbart_thr_tss_train, wbart_metrics_tss_train$tss, wbart_metrics_tss_train$fnr, wbart_metrics_tss_train$fpr))
  
  cat("\n=== dbart (Train) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_dbart_train, dbart_thr_fnr0_train, dbart_metrics_fnr0_train$tss, dbart_metrics_fnr0_train$fnr, dbart_metrics_fnr0_train$fpr,
              target_fnr, dbart_thr_fnr05_train, dbart_metrics_fnr05_train$tss, dbart_metrics_fnr05_train$fnr, dbart_metrics_fnr05_train$fpr,
              dbart_thr_tss_train, dbart_metrics_tss_train$tss, dbart_metrics_tss_train$fnr, dbart_metrics_tss_train$fpr))
  
  cat("\n=== wbart (Test) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_wbart_test, wbart_thr_fnr0, wbart_metrics_fnr0$tss, wbart_metrics_fnr0$fnr, wbart_metrics_fnr0$fpr,
              target_fnr, wbart_thr_fnr05, wbart_metrics_fnr05$tss, wbart_metrics_fnr05$fnr, wbart_metrics_fnr05$fpr,
              wbart_thr_tss, wbart_metrics_tss$tss, wbart_metrics_tss$fnr, wbart_metrics_tss$fpr))
  
  cat("\n=== dbart (Test) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold NTA = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_dbart_test, dbart_thr_fnr0, dbart_metrics_fnr0$tss, dbart_metrics_fnr0$fnr, dbart_metrics_fnr0$fpr,
              target_fnr, dbart_thr_fnr05, dbart_metrics_fnr05$tss, dbart_metrics_fnr05$fnr, dbart_metrics_fnr05$fpr,
              dbart_thr_tss, dbart_metrics_tss$tss, dbart_metrics_tss$fnr, dbart_metrics_tss$fpr))
  
  cat("\n=== bart.step (Test Quantile) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_bart_test_quant, bart_thr_fnr0_test_quant, bart_metrics_fnr0_test_quant$tss, bart_metrics_fnr0_test_quant$fnr, bart_metrics_fnr0_test_quant$fpr,
              target_fnr, bart_thr_fnr05_test_quant, bart_metrics_fnr05_test_quant$tss, bart_metrics_fnr05_test_quant$fnr, bart_metrics_fnr05_test_quant$fpr,
              bart_thr_tss_test_quant, bart_metrics_tss_test_quant$tss, bart_metrics_tss_test_quant$fnr, bart_metrics_tss_test_quant$fpr))
  
  cat("\n=== bart.step (Test GLM) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_bart_test_glm, bart_thr_fnr0_test_glm, bart_metrics_fnr0_test_glm$tss, bart_metrics_fnr0_test_glm$fnr, bart_metrics_fnr0_test_glm$fpr,
              target_fnr, bart_thr_fnr05_test_glm, bart_metrics_fnr05_test_glm$tss, bart_metrics_fnr05_test_glm$fnr, bart_metrics_fnr05_test_glm$fpr,
              bart_thr_tss_test_glm, bart_metrics_tss_test_glm$tss, bart_metrics_tss_test_glm$fnr, bart_metrics_tss_test_glm$fpr))
  
  cat("\n=== bart.step (Test LM) ===\n")
  cat(sprintf("AUC: %.4f\nFNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nFNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\nMax TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n",
              auc_bart_test_lm, bart_thr_fnr0_test_lm, bart_metrics_fnr0_test_lm$tss, bart_metrics_fnr0_test_lm$fnr, bart_metrics_fnr0_test_lm$fpr,
              target_fnr, bart_thr_fnr05_test_lm, bart_metrics_fnr05_test_lm$tss, bart_metrics_fnr05_test_lm$fnr, bart_metrics_fnr05_test_lm$fpr,
              bart_thr_tss_test_lm, bart_metrics_tss_test_lm$tss, bart_metrics_tss_test_lm$fnr, bart_metrics_tss_test_lm$fpr))
  
  # Variable Importance
  cat("\n=== Variable Importance (bart.step) ===\n")
  var_imp <- varimp(bart_model, plots = FALSE)
  if (is.null(var_imp)) {
    cat("Variable importance could not be computed.\n")
  } else {
    if (is.data.frame(var_imp) && "names" %in% colnames(var_imp) && "varimps" %in% colnames(var_imp)) {
      var_imp_df <- data.frame(
        Variable = as.character(var_imp$names),
        Importance = as.numeric(var_imp$varimps)
      )
      var_imp_df$Importance <- round(var_imp_df$Importance, 4)
      var_imp_df <- var_imp_df[order(-var_imp_df$Importance), ]
      cat(sprintf("%-30s %-10s\n", "Variable", "Importance"))
      cat(strrep("-", 40), "\n")
      for (i in 1:nrow(var_imp_df)) {
        cat(sprintf("%-30s %-10.4f\n", 
                    var_imp_df$Variable[i], 
                    var_imp_df$Importance[i]))
      }
      cat(strrep("-", 40), "\n")
    } else {
      cat("Variable importance data not found in varimp output.\n")
      cat("Structure of var_imp:\n")
      print(str(var_imp))
    }
  }
  
  # Determine Best Model Using AUC and KS-Test
  set.seed(123)
  bart_train_pred_jitter <- bart_train_pred + runif(length(bart_train_pred), -1e-6, 1e-6)
  bart_test_pred_quant_jitter <- bart_test_pred_quant + runif(length(bart_test_pred_quant), -1e-6, 1e-6)
  bart_test_pred_glm_jitter <- bart_test_pred_glm + runif(length(bart_test_pred_glm), -1e-6, 1e-6)
  bart_test_pred_lm_jitter <- bart_test_pred_lm + runif(length(bart_test_pred_lm), -1e-6, 1e-6)
  
  ks_quant <- ks.test(bart_train_pred_jitter, bart_test_pred_quant_jitter)$p.value
  ks_glm <- ks.test(bart_train_pred_jitter, bart_test_pred_glm_jitter)$p.value
  ks_lm <- ks.test(bart_train_pred_jitter, bart_test_pred_lm_jitter)$p.value
  
  auc_diff_quant <- abs(auc_bart_train - auc_bart_test_quant)
  auc_diff_glm <- abs(auc_bart_train - auc_bart_test_glm)
  auc_diff_lm <- abs(auc_bart_train - auc_bart_test_lm)
  
  model_metrics <- data.frame(
    Model = c("bart.step (Test Quantile)", "bart.step (Test GLM)", "bart.step (Test LM)"),
    KS_pvalue = c(ks_quant, ks_glm, ks_lm),
    AUC_diff = c(auc_diff_quant, auc_diff_glm, auc_diff_lm),
    AUC_test = c(auc_bart_test_quant, auc_bart_test_glm, auc_bart_test_lm)
  )
  
  max_auc_diff <- max(model_metrics$AUC_diff)
  max_ks_pvalue <- max(model_metrics$KS_pvalue)
  
  model_metrics$AUC_score <- if (max_auc_diff > 0) (1 - model_metrics$AUC_diff / max_auc_diff) else rep(1, nrow(model_metrics))
  model_metrics$KS_score <- if (max_ks_pvalue > 0) (model_metrics$KS_pvalue / max_ks_pvalue) else rep(1, nrow(model_metrics))
  
  model_metrics$Combined_score <- auc_weight * model_metrics$AUC_score + ks_weight * model_metrics$KS_score
  
  best_match_row <- model_metrics[which.max(model_metrics$Combined_score), ]
  best_match <- best_match_row$Model
  
  cat("\n=== Best Match for bart.step (Test) ===\n")
  cat("Best match (Weighted AUC + KS-test):", best_match, "\n")
  cat("Weights: AUC =", auc_weight, ", KS =", ks_weight, "\n")
  cat("KS p-values: Quantile =", ks_quant, ", GLM =", ks_glm, ", LM =", ks_lm, "\n")
  cat("AUC differences: Quantile =", auc_diff_quant, ", GLM =", auc_diff_glm, ", LM =", auc_diff_lm, "\n")
  cat("AUC (Test): Quantile =", auc_bart_test_quant, ", GLM =", auc_bart_test_glm, ", LM =", auc_bart_test_lm, "\n")
  cat("Combined scores:\n")
  cat(sprintf("%-30s %-10s\n", "Model", "Score"))
  cat(strrep("-" , 40), "\n")
  for (i in 1:nrow(model_metrics)) {
    cat(sprintf("%-30s %-10.4f\n", model_metrics$Model[i], model_metrics$Combined_score[i]))
  }
  cat(strrep("-", 40), "\n")
  
  # Select Best Test Model Metrics
  if (best_match == "bart.step (Test Quantile)") {
    best_test_pred <- bart_test_pred_quant
    best_test_auc <- auc_bart_test_quant
    best_test_metrics_fnr0 <- bart_metrics_fnr0_test_quant
    best_test_metrics_fnr05 <- bart_metrics_fnr05_test_quant
    best_test_metrics_tss <- bart_metrics_tss_test_quant
    best_test_thr_fnr0 <- bart_thr_fnr0_test_quant
    best_test_thr_fnr05 <- bart_thr_fnr05_test_quant
    best_test_thr_tss <- bart_thr_tss_test_quant
  } else if (best_match == "bart.step (Test GLM)") {
    best_test_pred <- bart_test_pred_glm
    best_test_auc <- auc_bart_test_glm
    best_test_metrics_fnr0 <- bart_metrics_fnr0_test_glm
    best_test_metrics_fnr05 <- bart_metrics_fnr05_test_glm
    best_test_metrics_tss <- bart_metrics_tss_test_glm
    best_test_thr_fnr0 <- bart_thr_fnr0_test_glm
    best_test_thr_fnr05 <- bart_thr_fnr05_test_glm
    best_test_thr_tss <- bart_thr_tss_test_glm
  } else {
    best_test_pred <- bart_test_pred_lm
    best_test_auc <- auc_bart_test_lm
    best_test_metrics_fnr0 <- bart_metrics_fnr0_test_lm
    best_test_metrics_fnr05 <- bart_metrics_fnr05_test_lm
    best_test_metrics_tss <- bart_metrics_tss_test_lm
    best_test_thr_fnr0 <- bart_thr_fnr0_test_lm
    best_test_thr_fnr05 <- bart_thr_fnr05_test_lm
    best_test_thr_tss <- bart_thr_tss_test_lm
  }
  
  # Overfitting Analysis
  # Overfitting Analysis
  cat("\n=== Overfitting Summary ===\n")
  auc_diff <- auc_bart_train - best_test_auc
  tss_diff <- bart_metrics_tss$tss - best_test_metrics_tss$tss
  fpr_diff <- bart_metrics_fnr05$fpr - best_test_metrics_fnr05$fpr
  ks_overfit <- ks.test(bart_train_pred, best_test_pred)
  residuals_train <- y.data - bart_train_pred
  residuals_test_best <- test_data$nf - best_test_pred
  ks_residuals <- ks.test(residuals_train, residuals_test_best)
  
  # Initialize overfitting_flags before spatial analysis
  overfitting_flags <- list()
  
  # Вычисление Moran's I для всех моделей
  if (do_spatial_analysis) {
    coords_train <- x.data[, c("x", "y")]
    coords_test <- test_data[, c("x", "y")]
    
    # Остатки для всех моделей
    residuals_bart_train <- y.data - bart_train_pred
    residuals_wbart_train <- y.data - wbart_train_pred
    residuals_dbart_train <- y.data - dbart_train_pred
    residuals_wbart_test <- test_data$nf - wbart_test_pred
    residuals_dbart_test <- test_data$nf - dbart_test_pred
    residuals_bart_test_quant <- test_data$nf - bart_test_pred_quant
    residuals_bart_test_glm <- test_data$nf - bart_test_pred_glm
    residuals_bart_test_lm <- test_data$nf - bart_test_pred_lm
    
    # Вычисление Moran's I для каждой модели
    optimal_k_bart <- find_optimal_k(residuals_bart_train, coords_train)
    moran_bart_train <- moran.test(residuals_bart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_bart$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_wbart_train <- find_optimal_k(residuals_wbart_train, coords_train)
    moran_wbart_train <- moran.test(residuals_wbart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_wbart_train$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_dbart_train <- find_optimal_k(residuals_dbart_train, coords_train)
    moran_dbart_train <- moran.test(residuals_dbart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_dbart_train$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_wbart_test <- find_optimal_k(residuals_wbart_test, coords_test)
    moran_wbart_test <- moran.test(residuals_wbart_test, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_wbart_test$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_dbart_test <- find_optimal_k(residuals_dbart_test, coords_test)
    moran_dbart_test <- moran.test(residuals_dbart_test, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_dbart_test$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_quant <- find_optimal_k(residuals_bart_test_quant, coords_test)
    moran_bart_test_quant <- moran.test(residuals_bart_test_quant, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_quant$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_glm <- find_optimal_k(residuals_bart_test_glm, coords_test)
    moran_bart_test_glm <- moran.test(residuals_bart_test_glm, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_glm$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_lm <- find_optimal_k(residuals_bart_test_lm, coords_test)
    moran_bart_test_lm <- moran.test(residuals_bart_test_lm, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_lm$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    cross_moran_result <- cross_moran(residuals_bart_train, residuals_test_best, coords_train, coords_test, k = optimal_k_bart_test_lm$optimal_k)
    
    moran_train <- moran_bart_train$estimate[1]
    moran_test <- moran_bart_test_lm$estimate[1]
    cross_moran <- cross_moran_result$estimate[1]
    
    # Флаги переобучения
    if (moran_train <= 0.2) {
      overfitting_flags$moran_train <- "Good: Train residuals Moran's I ≤ 0.2, model accounts for spatial autocorrelation"
    } else if (moran_train > 0.3) {
      overfitting_flags$moran_train <- "Suspicious: Train residuals Moran's I > 0.3, model fails to account for spatial autocorrelation"
    } else {
      overfitting_flags$moran_train <- "Moderate: Train residuals Moran's I between 0.2 and 0.3, partial accounting for spatial autocorrelation"
    }
    
    if (moran_test <= 0.2) {
      overfitting_flags$moran_test <- "Good: Test residuals Moran's I ≤ 0.2, no overfitting to spatial noise"
    } else if (moran_test > 0.3) {
      overfitting_flags$moran_test <- "Suspicious: Test residuals Moran's I > 0.3, potential overfitting to spatial noise"
    } else {
      overfitting_flags$moran_test <- "Moderate: Test residuals Moran's I between 0.2 and 0.3, possible minor overfitting to spatial patterns"
    }
    
    if (cross_moran <= 0.2) {
      overfitting_flags$cross_moran <- "Good: Cross-Moran's I ≤ 0.2, spatial patterns independent between train and test"
    } else {
      overfitting_flags$cross_moran <- "Suspicious: Cross-Moran's I > 0.2, potential overfitting to shared spatial patterns"
    }
  }
  
  # Вывод обновленной таблицы Moran's I
  # Overfitting Analysis
  cat("\n=== Overfitting Summary ===\n")
  auc_diff <- auc_bart_train - best_test_auc
  tss_diff <- bart_metrics_tss$tss - best_test_metrics_tss$tss
  fpr_diff <- bart_metrics_fnr05$fpr - best_test_metrics_fnr05$fpr
  ks_overfit <- ks.test(bart_train_pred, best_test_pred)
  residuals_train <- y.data - bart_train_pred
  residuals_test_best <- test_data$nf - best_test_pred
  ks_residuals <- ks.test(residuals_train, residuals_test_best)
  
  # Вычисление Moran's I для всех моделей с оптимальным k
  if (do_spatial_analysis) {
    coords_train <- x.data[, c("x", "y")]
    coords_test <- test_data[, c("x", "y")]
    
    # Остатки для всех моделей
    residuals_bart_train <- y.data - bart_train_pred
    residuals_wbart_train <- y.data - wbart_train_pred
    residuals_dbart_train <- y.data - dbart_train_pred
    residuals_wbart_test <- test_data$nf - wbart_test_pred
    residuals_dbart_test <- test_data$nf - dbart_test_pred
    residuals_bart_test_quant <- test_data$nf - bart_test_pred_quant
    residuals_bart_test_glm <- test_data$nf - bart_test_pred_glm
    residuals_bart_test_lm <- test_data$nf - bart_test_pred_lm
    
    # Вычисление оптимального k и Moran's I для каждой модели
    optimal_k_bart <- find_optimal_k(residuals_bart_train, coords_train)
    moran_bart_train <- moran.test(residuals_bart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_bart$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_wbart_train <- find_optimal_k(residuals_wbart_train, coords_train)
    moran_wbart_train <- moran.test(residuals_wbart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_wbart_train$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_dbart_train <- find_optimal_k(residuals_dbart_train, coords_train)
    moran_dbart_train <- moran.test(residuals_dbart_train, nb2listw(knn2nb(knearneigh(coords_train, k = optimal_k_dbart_train$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_wbart_test <- find_optimal_k(residuals_wbart_test, coords_test)
    moran_wbart_test <- moran.test(residuals_wbart_test, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_wbart_test$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_dbart_test <- find_optimal_k(residuals_dbart_test, coords_test)
    moran_dbart_test <- moran.test(residuals_dbart_test, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_dbart_test$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_quant <- find_optimal_k(residuals_bart_test_quant, coords_test)
    moran_bart_test_quant <- moran.test(residuals_bart_test_quant, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_quant$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_glm <- find_optimal_k(residuals_bart_test_glm, coords_test)
    moran_bart_test_glm <- moran.test(residuals_bart_test_glm, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_glm$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    optimal_k_bart_test_lm <- find_optimal_k(residuals_bart_test_lm, coords_test)
    moran_bart_test_lm <- moran.test(residuals_bart_test_lm, nb2listw(knn2nb(knearneigh(coords_test, k = optimal_k_bart_test_lm$optimal_k)), style = "W", zero.policy = TRUE), zero.policy = TRUE)
    
    # Кросс-корреляция с использованием k от лучшей тестовой модели (bart.step (Test LM))
    cross_moran_result <- cross_moran(residuals_bart_train, residuals_test_best, coords_train, coords_test, k = optimal_k_bart_test_lm$optimal_k)
    
    # Формирование таблицы Moran's I
    cat("\n=== Summary Table of Moran's I ===\n")
    cat(sprintf("%-30s %-10s %-10s %-10s\n", "Model", "k", "Moran I", "P-Value"))
    cat(strrep("-", 60), "\n")
    moran_summary <- data.frame(
      Model = c("bart.step (Train)", "wbart (Train)", "dbart (Train)", 
                "wbart (Test)", "dbart (Test)", 
                "bart.step (Test Quantile)", "bart.step (Test GLM)", "bart.step (Test LM)", 
                "bart.step (Train vs Test)"),
      k = c(optimal_k_bart$optimal_k, optimal_k_wbart_train$optimal_k, optimal_k_dbart_train$optimal_k,
            optimal_k_wbart_test$optimal_k, optimal_k_dbart_test$optimal_k,
            optimal_k_bart_test_quant$optimal_k, optimal_k_bart_test_glm$optimal_k, optimal_k_bart_test_lm$optimal_k,
            optimal_k_bart_test_lm$optimal_k),
      Moran_I = c(moran_bart_train$estimate[1], moran_wbart_train$estimate[1], moran_dbart_train$estimate[1],
                  moran_wbart_test$estimate[1], moran_dbart_test$estimate[1],
                  moran_bart_test_quant$estimate[1], moran_bart_test_glm$estimate[1], moran_bart_test_lm$estimate[1],
                  cross_moran_result$estimate[1]),
      P_Value = c(moran_bart_train$p.value, moran_wbart_train$p.value, moran_dbart_train$p.value,
                  moran_wbart_test$p.value, moran_dbart_test$p.value,
                  moran_bart_test_quant$p.value, moran_bart_test_glm$p.value, moran_bart_test_lm$p.value,
                  cross_moran_result$p.value)
    )
    for (i in 1:nrow(moran_summary)) {
      p_val <- moran_summary$P_Value[i]
      p_val_display <- if (p_val < 0.0001) "< 1e-04" else sprintf("%.6f", p_val)
      cat(sprintf("%-30s %-10d %-10.4f %-10s\n", 
                  moran_summary$Model[i], 
                  moran_summary$k[i], 
                  moran_summary$Moran_I[i], 
                  p_val_display))
    }
    cat(strrep("-", 60), "\n")
  } else {
    cat("\nSpatial analysis not performed, so Moran's I summary is unavailable.\n")
  }
  
  # Формирование флагов переобучения с использованием данных из moran_summary
  overfitting_flags <- list()
  
  # AUC difference
  if (auc_diff < 0.05) {
    overfitting_flags$auc <- sprintf("AUC difference (%.4f): Normal: AUC difference < 0.05", auc_diff)
  } else if (auc_diff >= 0.05 && auc_diff <= 0.10) {
    overfitting_flags$auc <- sprintf("AUC difference (%.4f): Possible overfitting: AUC difference between 0.05 and 0.10", auc_diff)
  } else {
    overfitting_flags$auc <- sprintf("AUC difference (%.4f): Likely overfitting: AUC difference > 0.10", auc_diff)
  }
  
  # TSS difference
  if (tss_diff < 0.1) {
    overfitting_flags$tss <- sprintf("TSS difference (%.4f): Normal: TSS difference < 0.1", tss_diff)
  } else if (tss_diff >= 0.1 && tss_diff <= 0.2) {
    overfitting_flags$tss <- sprintf("TSS difference (%.4f): Suspicious: TSS difference between 0.1 and 0.2", tss_diff)
  } else {
    overfitting_flags$tss <- sprintf("TSS difference (%.4f): Overfitting: TSS difference > 0.2", tss_diff)
  }
  
  # FPR difference
  if (abs(fpr_diff) < 0.05) {
    overfitting_flags$fpr <- sprintf("FPR difference at FNR=0.0400 (%.4f): Normal: FPR difference within acceptable range", fpr_diff)
  } else {
    overfitting_flags$fpr <- sprintf("FPR difference at FNR=0.0400 (%.4f): Suspicious: FPR difference outside acceptable range", fpr_diff)
  }
  
  # KS-Test for Predicted Probabilities
  if (ks_overfit$p.value >= 0.05) {
    overfitting_flags$ks_prob <- sprintf("KS-Test for Predicted Probabilities (p-value = %.6f): Good generalization: No significant difference in predicted probability distributions (KS p-value >= 0.05)", ks_overfit$p.value)
  } else {
    overfitting_flags$ks_prob <- sprintf("KS-Test for Predicted Probabilities (p-value = %.6f): Potential overfitting: Significant difference in predicted probability distributions (KS p-value < 0.05)", ks_overfit$p.value)
  }
  
  # KS-Test for Residuals
  if (ks_residuals$p.value >= 0.05) {
    overfitting_flags$ks_resid <- sprintf("KS-Test for Residuals (p-value = %.6f): Good generalization: No significant difference in residual distributions (KS p-value >= 0.05)", ks_residuals$p.value)
  } else {
    overfitting_flags$ks_resid <- sprintf("KS-Test for Residuals (p-value = %.6f): Potential overfitting: Significant difference in residual distributions (KS p-value < 0.05)", ks_residuals$p.value)
  }
  
  # Moran's I flags using moran_summary
  if (do_spatial_analysis) {
    # Moran's I for Train (bart.step (Train))
    moran_train <- moran_summary$Moran_I[moran_summary$Model == "bart.step (Train)"]
    moran_train_pval <- moran_summary$P_Value[moran_summary$Model == "bart.step (Train)"]
    moran_train_k <- moran_summary$k[moran_summary$Model == "bart.step (Train)"]
    if (moran_train <= 0.2) {
      overfitting_flags$moran_train <- sprintf("Moran's I for Train Residuals (%.4f, k=%d, p-value=%.6f): Good: Train residuals Moran's I ≤ 0.2, model accounts for spatial autocorrelation", moran_train, moran_train_k, moran_train_pval)
    } else if (moran_train > 0.3) {
      overfitting_flags$moran_train <- sprintf("Moran's I for Train Residuals (%.4f, k=%d, p-value=%.6f): Suspicious: Train residuals Moran's I > 0.3, model fails to account for spatial autocorrelation", moran_train, moran_train_k, moran_train_pval)
    } else {
      overfitting_flags$moran_train <- sprintf("Moran's I for Train Residuals (%.4f, k=%d, p-value=%.6f): Moderate: Train residuals Moran's I between 0.2 and 0.3, partial accounting for spatial autocorrelation", moran_train, moran_train_k, moran_train_pval)
    }
    
    # Moran's I for Test (best model, e.g., bart.step (Test LM))
    moran_test <- moran_summary$Moran_I[moran_summary$Model == "bart.step (Test LM)"]
    moran_test_pval <- moran_summary$P_Value[moran_summary$Model == "bart.step (Test LM)"]
    moran_test_k <- moran_summary$k[moran_summary$Model == "bart.step (Test LM)"]
    if (moran_test <= 0.2) {
      overfitting_flags$moran_test <- sprintf("Moran's I for Test Residuals (%.4f, k=%d, p-value=%.6f): Good: Test residuals Moran's I ≤ 0.2, no overfitting to spatial noise", moran_test, moran_test_k, moran_test_pval)
    } else if (moran_test > 0.3) {
      overfitting_flags$moran_test <- sprintf("Moran's I for Test Residuals (%.4f, k=%d, p-value=%.6f): Suspicious: Test residuals Moran's I > 0.3, potential overfitting to spatial noise", moran_test, moran_test_k, moran_test_pval)
    } else {
      overfitting_flags$moran_test <- sprintf("Moran's I for Test Residuals (%.4f, k=%d, p-value=%.6f): Moderate: Test residuals Moran's I between 0.2 and 0.3, possible minor overfitting to spatial patterns", moran_test, moran_test_k, moran_test_pval)
    }
    
    # Cross-Moran's I
    cross_moran <- moran_summary$Moran_I[moran_summary$Model == "bart.step (Train vs Test)"]
    cross_moran_pval <- moran_summary$P_Value[moran_summary$Model == "bart.step (Train vs Test)"]
    cross_moran_k <- moran_summary$k[moran_summary$Model == "bart.step (Train vs Test)"]
    if (cross_moran <= 0.2) {
      overfitting_flags$cross_moran <- sprintf("Cross-Moran's I (%.4f, k=%d, p-value=%.6f): Good: Cross-Moran's I ≤ 0.2, spatial patterns independent between train and test", cross_moran, cross_moran_k, cross_moran_pval)
    } else {
      overfitting_flags$cross_moran <- sprintf("Cross-Moran's I (%.4f, k=%d, p-value=%.6f): Suspicious: Cross-Moran's I > 0.2, potential overfitting to shared spatial patterns", cross_moran, cross_moran_k, cross_moran_pval)
    }
  }
  
  # Вывод результатов анализа переобучения
  cat("\nOverfitting Analysis Results:\n")
  for (flag in overfitting_flags) {
    cat(sprintf("- %s\n", flag))
  }
  
  # Вывод заключения о переобучении
  cat("\nConclusion: ")
  if (any(grepl("Suspicious|Potential overfitting|Likely overfitting|Overfitting", overfitting_flags))) {
    cat("Potential overfitting detected based on the following:\n")
    for (flag in overfitting_flags) {
      if (grepl("Suspicious|Potential overfitting|Likely overfitting|Overfitting", flag)) {
        cat(sprintf("- %s\n", flag))
      }
    }
    if (do_spatial_analysis && (moran_train > 0.3 || moran_test > 0.3 || cross_moran > 0.2)) {
      cat("- Spatial overfitting likely due to variables like wc2.1_10m_elev capturing spatial patterns.\n")
    }
  } else {
    cat("No significant overfitting detected.\n") 
  }
  
  # Сохранение гистограмм остатков
  if (do_spatial_analysis) {
    pdf("residuals_histograms.pdf")
    hist(residuals_train, main = "Histogram of Train Residuals", xlab = "Residuals", col = "lightblue")
    hist(residuals_test_best, main = "Histogram of Test Residuals (Best Model)", xlab = "Residuals", col = "lightgreen")
    dev.off()
    cat("Residual histograms saved to residuals_histograms.pdf\n")
  }
  # Output Overfitting Analysis
  cat("Overfitting Analysis Results:\n")
  cat(sprintf("- AUC difference (%.4f): %s\n", auc_diff, overfitting_flags$auc))
  cat(sprintf("- TSS difference (%.4f): %s\n", tss_diff, overfitting_flags$tss))
  cat(sprintf("- FPR difference at FNR=%.4f (%.4f): %s\n", target_fnr, fpr_diff, 
              ifelse(abs(fpr_diff) > 0.2, "Suspicious: Large FPR difference (> 0.2)", "Normal: FPR difference within acceptable range")))
  cat(sprintf("- KS-Test for Predicted Probabilities (p-value = %.6f): %s\n", ks_overfit$p.value, overfitting_flags$ks_prob))
  cat(sprintf("- KS-Test for Residuals (p-value = %.6f): %s\n", ks_residuals$p.value, overfitting_flags$ks_resid))
  
  if (do_spatial_analysis) {
    cat(sprintf("- Moran's I for Train Residuals (%.4f, k=%d, p-value=%.6f): %s\n", 
                moran_train, optimal_k_bart$optimal_k, moran_bart_train$p.value, overfitting_flags$moran_train))
    cat(sprintf("- Moran's I for Test Residuals (%.4f, k=%d, p-value=%.6f): %s\n", 
                moran_test, optimal_k_bart_test_lm$optimal_k, moran_bart_test_lm$p.value, overfitting_flags$moran_test))
    cat(sprintf("- Cross-Moran's I (%.4f, k=%d, p-value=%.6f): %s\n", 
                cross_moran, optimal_k_bart_test_lm$optimal_k, cross_moran_result$p.value, overfitting_flags$cross_moran))
  }
  
  overfitting_detected <- any(grepl("overfitting|Suspicious", unlist(overfitting_flags)))
  if (overfitting_detected) {
    cat("\nConclusion: Potential overfitting detected based on the following:\n")
    for (flag_name in names(overfitting_flags)) {
      if (grepl("overfitting|Suspicious", overfitting_flags[[flag_name]])) {
        cat(sprintf("- %s\n", overfitting_flags[[flag_name]]))
      }
    }
    if (do_spatial_analysis && (moran_test > 0.3 || cross_moran > 0.2)) {
      cat("- Spatial overfitting likely due to variables like wc2.1_10m_elev capturing spatial patterns.\n")
    }
  } else {
    cat("\nConclusion: No strong evidence of overfitting. All metrics within acceptable ranges.\n")
  }
  
  # Residual Histograms
  residuals_data <- data.frame(
    Residuals = c(residuals_train, residuals_test_best),
    Set = factor(rep(c("bart.step (Train)", best_match), 
                     times = c(length(residuals_train), length(residuals_test_best)))),
    True_Class = factor(c(y.data, test_data$nf), levels = c(0, 1), 
                        labels = c("Negative (0)", "Positive (1)"))
  )
  
  p_residuals <- ggplot(residuals_data, aes(x = Residuals, fill = True_Class)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    facet_wrap(~Set, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Negative (0)" = "gray", "Positive (1)" = "orange")) +
    labs(title = "Distribution of Residuals (Train vs Best Test Model)",
         x = "Residuals", y = "Count", fill = "True Class") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  pdf("residuals_histograms.pdf", width = 8, height = 6)
  print(p_residuals)
  dev.off()
  cat("Residual histograms saved to residuals_histograms.pdf\n")
  
  # Prediction Histograms
  pred_data <- data.frame(
    Probability = c(bart_train_pred, wbart_train_pred, dbart_train_pred, wbart_test_pred, dbart_test_pred, bart_test_pred_quant, bart_test_pred_glm, bart_test_pred_lm),
    True_Class = c(y.data, y.data, y.data, test_data$nf, test_data$nf, test_data$nf, test_data$nf, test_data$nf),
    Set = factor(rep(c("bart.step (Train)", "wbart (Train)", "dbart (Train)", "wbart (Test)", "dbart (Test)", "bart.step (Test Quantile)", "bart.step (Test GLM)", "bart.step (Test LM)"), 
                     times = c(length(bart_train_pred), length(wbart_train_pred), length(dbart_train_pred), length(wbart_test_pred), length(dbart_test_pred), length(bart_test_pred_quant), length(bart_test_pred_glm), length(bart_test_pred_lm))),
                 levels = c(best_match, setdiff(c("bart.step (Train)", "wbart (Train)", "dbart (Train)", "wbart (Test)", "dbart (Test)", "bart.step (Test Quantile)", "bart.step (Test GLM)", "bart.step (Test LM)"), best_match)))
  )
  
  pred_data$Probability <- pmin(pmax(pred_data$Probability, 0), 1)
  pred_data$True_Class <- factor(pred_data$True_Class, levels = c(0, 1), labels = c("Negative (0)", "Positive (1)"))
  
  cat("\nChecking pred_data for NA before plotting:\n")
  cat("NA in Probability:", sum(is.na(pred_data$Probability)), "\n")
  cat("NA in True_Class:", sum(is.na(pred_data$True_Class)), "\n")
  pred_data <- pred_data[!is.na(pred_data$Probability) & !is.na(pred_data$True_Class), ]
  
  out_of_range <- sum(pred_data$Probability < 0 | pred_data$Probability > 1)
  cat("Values outside [0, 1] range:", out_of_range, "\n")
  
  thresholds <- data.frame(
    Set = factor(c("bart.step (Train)", "wbart (Train)", "dbart (Train)", "wbart (Test)", "dbart (Test)", "bart.step (Test Quantile)", "bart.step (Test GLM)", "bart.step (Test LM)"),
                 levels = levels(pred_data$Set)),
    FNR0 = c(bart_thr_fnr0, wbart_thr_fnr0_train, dbart_thr_fnr0_train, wbart_thr_fnr0, dbart_thr_fnr0, bart_thr_fnr0_test_quant, bart_thr_fnr0_test_glm, bart_thr_fnr0_test_lm),
    FNR05 = c(bart_thr_fnr05, wbart_thr_fnr05_train, dbart_thr_fnr05_train, wbart_thr_fnr05, dbart_thr_fnr05, bart_thr_fnr05_test_quant, bart_thr_fnr05_test_glm, bart_thr_fnr05_test_lm),
    TSS = c(bart_thr_tss, wbart_thr_tss_train, dbart_thr_tss_train, wbart_thr_tss, dbart_thr_tss, bart_thr_tss_test_quant, bart_thr_tss_test_glm, bart_thr_tss_test_lm)
  )
  
  fill_data <- data.frame(
    Set = thresholds$Set,
    xmin = thresholds$FNR0,
    xmax = thresholds$FNR05,
    ymin = -Inf,
    ymax = Inf
  )
  
  thresholds_long <- reshape2::melt(thresholds, id.vars = "Set", measure.vars = c("FNR0", "FNR05", "TSS"), variable.name = "Threshold_Type", value.name = "Value")
  thresholds_long$Threshold_Type <- factor(thresholds_long$Threshold_Type, levels = c("FNR0", "FNR05", "TSS"), labels = c("FNR=0", "FNR=0.05", "TSS"))
  
  cat("Checking thresholds_long for NA:\n")
  cat("NA in Value:", sum(is.na(thresholds_long$Value)), "\n")
  thresholds_long <- thresholds_long[!is.na(thresholds_long$Value), ]
  # Проверка типов данных в thresholds_long
  cat("Проверка типов данных в thresholds_long:\n")
  cat("Тип данных Value:", class(thresholds_long$Value), "\n")
  cat("Тип данных Threshold_Type:", class(thresholds_long$Threshold_Type), "\n")
  
  # Проверка значений linewidth
  linewidth_values <- c("FNR=0" = 0.5, "FNR=0.05" = 0.5, "TSS" = 1.5)
  if (!all(sapply(linewidth_values, is.numeric))) {
    stop("Все значения в scale_linewidth_manual должны быть числовыми!")
  }
  # Создание гистограммы предсказанных вероятностей
  p_hist <- ggplot(pred_data, aes(x = Probability, fill = True_Class)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity", boundary = 0) +
    facet_wrap(~Set, scales = "free_y", ncol = 2, strip.position = "top", drop = FALSE) +
    geom_rect(data = fill_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
    geom_vline(data = thresholds_long, aes(xintercept = Value, color = Threshold_Type, linetype = Threshold_Type, linewidth = Threshold_Type), 
               show.legend = TRUE) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.05),
      limits = c(0, 1),
      name = "Probability",
      sec.axis = sec_axis(~ ., breaks = seq(0, 1, by = 0.05), name = "Probability (Top)")
    ) +
    scale_fill_manual(values = c("Negative (0)" = "gray", "Positive (1)" = "orange")) +
    scale_color_manual(values = c("FNR=0" = "blue", "FNR=0.05" = "red", "TSS" = "green")) +
    scale_linetype_manual(values = c("FNR=0" = "solid", "FNR=0.05" = "dashed", "TSS" = "dotted")) +
    scale_linewidth_manual(values = c("FNR=0" = 0.5, "FNR=0.05" = 0.5, "TSS" = 1.5)) +  # Исправлено: все значения числовые
    labs(title = paste("Distribution of Predicted Probabilities (Best Model:", best_match, ")"), 
         x = "Probability", y = "Count", 
         fill = "True Class", color = "Threshold", linetype = "Threshold", linewidth = "Threshold") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.title.x.top = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  pdf("prediction_histograms.pdf", width = 12, height = 10)
  print(p_hist)
  dev.off()
  cat("Histograms saved to prediction_histograms.pdf\n")
  
  # Final Summary
  cat("\n=== Final Summary ===\n")
  
  cat("\n=== ", best_match, " ===\n")
  cat(sprintf("AUC: %.4f\n", best_test_auc))
  cat(sprintf("FNR=0 (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n", 
              best_test_thr_fnr0, best_test_metrics_fnr0$tss, best_test_metrics_fnr0$fnr, best_test_metrics_fnr0$fpr))
  cat(sprintf("FNR=%.4f (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n", 
              target_fnr, best_test_thr_fnr05, best_test_metrics_fnr05$tss, best_test_metrics_fnr05$fnr, best_test_metrics_fnr05$fpr))
  cat(sprintf("Max TSS (threshold = %.4f): TSS = %.4f, FNR = %.4f, FPR = %.4f\n", 
              best_test_thr_tss, best_test_metrics_tss$tss, best_test_metrics_tss$fnr, best_test_metrics_tss$fpr))
  
  cat("Predictor list:", paste(final_vars, collapse = " "), "\n")
  
  cat("\n=== Overfitting Summary ===\n")
  if (overfitting_detected) {
    cat("Potential overfitting detected:\n")
    for (flag_name in names(overfitting_flags)) {
      if (grepl("overfitting|Suspicious", overfitting_flags[[flag_name]])) {
        cat(sprintf("- %s\n", overfitting_flags[[flag_name]]))
      }
    }
    if (do_spatial_analysis && (moran_test > 0.3 || cross_moran > 0.2)) {
      cat("- Spatial overfitting likely due to variables like wc2.1_10m_elev capturing spatial patterns.\n")
    }
  } else {
    cat("No strong evidence of overfitting:\n")
    cat(sprintf("- AUC difference (%.4f): %s\n", auc_diff, overfitting_flags$auc))
    cat(sprintf("- TSS difference (%.4f): %s\n", tss_diff, overfitting_flags$tss))
    cat(sprintf("- KS-Test for Predicted Probabilities (p-value = %.6f): %s\n", ks_overfit$p.value, overfitting_flags$ks_prob))
    cat(sprintf("- KS-Test for Residuals (p-value = %.6f): %s\n", ks_residuals$p.value, overfitting_flags$ks_resid))
  }
  
  if (do_spatial_analysis) {
    cat("\n=== Summary Table of Moran's I ===\n")
    cat(sprintf("%-30s %-10s %-10s %-10s\n", "Model", "k", "Moran I", "P-Value"))
    cat(strrep("-", 60), "\n")
    moran_summary <- data.frame(
      Model = c("bart.step (Train)", best_match, "bart.step (Train vs Test)"),
      k = c(optimal_k_bart$optimal_k, optimal_k_bart_test_lm$optimal_k, optimal_k_bart_test_lm$optimal_k),
      Moran_I = c(moran_bart_train$estimate[1], moran_bart_test_lm$estimate[1], cross_moran_result$estimate[1]),
      P_Value = c(moran_bart_train$p.value, moran_bart_test_lm$p.value, cross_moran_result$p.value)
    )
    for (i in 1:nrow(moran_summary)) {
      p_val <- moran_summary$P_Value[i]
      p_val_display <- if (p_val < 0.0001) "< 1e-04" else sprintf("%.6f", p_val)
      cat(sprintf("%-30s %-10d %-10.4f %-10s\n", 
                  moran_summary$Model[i], 
                  moran_summary$k[i], 
                  moran_summary$Moran_I[i], 
                  p_val_display))
    }
    cat(strrep("-", 60), "\n")
  } else {
    cat("\nSpatial analysis not performed, so Moran's I summary is unavailable.\n")
  }
  
  # Statistical Summary Table for Variables
  cat("\n=== Statistical Summary of Variables ===\n")
  
  # Split data into groups based on y.data (0 and 1)
  group0 <- x.data[y.data == 0, final_vars, drop = FALSE]
  group1 <- x.data[y.data == 1, final_vars, drop = FALSE]
  
  # Calculate statistics for each variable
  stats_summary <- data.frame(
    Variable = final_vars,
    Importance = NA,
    M_Group0 = NA,
    SD_Group0 = NA,
    SE_Group0 = NA,
    M_Group1 = NA,
    SD_Group1 = NA,
    SE_Group1 = NA,
    T_Statistic = NA,
    P_Value = NA
  )
  
  # Get variable importance
  if (!is.null(var_imp) && is.data.frame(var_imp) && "names" %in% colnames(var_imp) && "varimps" %in% colnames(var_imp)) {
    for (i in 1:nrow(stats_summary)) {
      var_name <- stats_summary$Variable[i]
      imp_idx <- which(var_imp$names == var_name)
      if (length(imp_idx) > 0) {
        stats_summary$Importance[i] <- var_imp$varimps[imp_idx]
      }
    }
  }
  
  # Calculate means, SD, SE, and t-test for each variable
  for (i in 1:nrow(stats_summary)) {
    var_name <- stats_summary$Variable[i]
    
    # Group 0 (y=0)
    data0 <- group0[[var_name]]
    stats_summary$M_Group0[i] <- mean(data0, na.rm = TRUE)
    stats_summary$SD_Group0[i] <- sd(data0, na.rm = TRUE)
    stats_summary$SE_Group0[i] <- stats_summary$SD_Group0[i] / sqrt(sum(!is.na(data0)))
    
    # Group 1 (y=1)
    data1 <- group1[[var_name]]
    stats_summary$M_Group1[i] <- mean(data1, na.rm = TRUE)
    stats_summary$SD_Group1[i] <- sd(data1, na.rm = TRUE)
    stats_summary$SE_Group1[i] <- stats_summary$SD_Group1[i] / sqrt(sum(!is.na(data1)))
    
    # t-test
    t_test <- tryCatch({
      t.test(data0, data1, var.equal = FALSE)
    }, error = function(e) {
      list(statistic = NA, p.value = NA)
    })
    stats_summary$T_Statistic[i] <- t_test$statistic
    stats_summary$P_Value[i] <- t_test$p.value
  }
  
  # Round numerical columns
  stats_summary$Importance <- round(stats_summary$Importance, 4)
  stats_summary$M_Group0 <- round(stats_summary$M_Group0, 3)
  stats_summary$SD_Group0 <- round(stats_summary$SD_Group0, 3)
  stats_summary$SE_Group0 <- round(stats_summary$SE_Group0, 3)
  stats_summary$M_Group1 <- round(stats_summary$M_Group1, 3)
  stats_summary$SD_Group1 <- round(stats_summary$SD_Group1, 3)
  stats_summary$SE_Group1 <- round(stats_summary$SE_Group1, 3)
  stats_summary$T_Statistic <- round(stats_summary$T_Statistic, 3)
  
  # Print the table
  cat(sprintf("%-20s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n",
              "Variable", "Importance", "M (y=0)", "SD (y=0)", "SE (y=0)", "M (y=1)", "SD (y=1)", "SE (y=1)", "t-stat", "p-value"))
  cat(strrep("-", 100), "\n")
  for (i in 1:nrow(stats_summary)) {
    p_val <- stats_summary$P_Value[i]
    p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.0001) "< 1e-04" else sprintf("%.4f", p_val)
    cat(sprintf("%-20s %-10.4f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10s\n",
                stats_summary$Variable[i],
                stats_summary$Importance[i],
                stats_summary$M_Group0[i],
                stats_summary$SD_Group0[i],
                stats_summary$SE_Group0[i],
                stats_summary$M_Group1[i],
                stats_summary$SD_Group1[i],
                stats_summary$SE_Group1[i],
                stats_summary$T_Statistic[i],
                p_val_display))
  }
  cat(strrep("-", 100), "\n")
  cat("Note: Importance is from bart.step variable importance; M, SD, SE, t-stat, and p-value are calculated for training data split by y.data (0 and 1).\n")
  
  # Return Results
  results <- list(
    bart_model = bart_model,
    wbart_model = wbart_model,
    dbart_model = dbart_model,
    final_vars = final_vars,
    bart_train_pred = bart_train_pred,
    wbart_train_pred = wbart_train_pred,
    dbart_train_pred = dbart_train_pred,
    wbart_test_pred = wbart_test_pred,
    dbart_test_pred = dbart_test_pred,
    bart_test_pred_quant = bart_test_pred_quant,
    bart_test_pred_glm = bart_test_pred_glm,
    bart_test_pred_lm = bart_test_pred_lm,
    auc = list(
      bart_train = auc_bart_train,
      wbart_train = auc_wbart_train,
      dbart_train = auc_dbart_train,
      wbart_test = auc_wbart_test,
      dbart_test = auc_dbart_test,
      bart_test_quant = auc_bart_test_quant,
      bart_test_glm = auc_bart_test_glm,
      bart_test_lm = auc_bart_test_lm
    ),
    metrics = list(
      bart_train = list(fnr0 = bart_metrics_fnr0, fnr05 = bart_metrics_fnr05, tss = bart_metrics_tss),
      wbart_train = list(fnr0 = wbart_metrics_fnr0_train, fnr05 = wbart_metrics_fnr05_train, tss = wbart_metrics_tss_train),
      dbart_train = list(fnr0 = dbart_metrics_fnr0_train, fnr05 = dbart_metrics_fnr05_train, tss = dbart_metrics_tss_train),
      wbart_test = list(fnr0 = wbart_metrics_fnr0, fnr05 = wbart_metrics_fnr05, tss = wbart_metrics_tss),
      dbart_test = list(fnr0 = dbart_metrics_fnr0, fnr05 = dbart_metrics_fnr05, tss = dbart_metrics_tss),
      bart_test_quant = list(fnr0 = bart_metrics_fnr0_test_quant, fnr05 = bart_metrics_fnr05_test_quant, tss = bart_metrics_tss_test_quant),
      bart_test_glm = list(fnr0 = bart_metrics_fnr0_test_glm, fnr05 = bart_metrics_fnr05_test_glm, tss = bart_metrics_tss_test_glm),
      bart_test_lm = list(fnr0 = bart_metrics_fnr0_test_lm, fnr05 = bart_metrics_fnr05_test_lm, tss = bart_metrics_tss_test_lm)
    ),
    thresholds = list(
      bart_train = list(fnr0 = bart_thr_fnr0, fnr05 = bart_thr_fnr05, tss = bart_thr_tss),
      wbart_train = list(fnr0 = wbart_thr_fnr0_train, fnr05 = wbart_thr_fnr05_train, tss = wbart_thr_tss_train),
      dbart_train = list(fnr0 = dbart_thr_fnr0_train, fnr05 = dbart_thr_fnr05_train, tss = dbart_thr_tss_train),
      wbart_test = list(fnr0 = wbart_thr_fnr0, fnr05 = wbart_thr_fnr05, tss = wbart_thr_tss),
      dbart_test = list(fnr0 = dbart_thr_fnr0, fnr05 = dbart_thr_fnr05, tss = dbart_thr_tss),
      bart_test_quant = list(fnr0 = bart_thr_fnr0_test_quant, fnr05 = bart_thr_fnr05_test_quant, tss = bart_thr_tss_test_quant),
      bart_test_glm = list(fnr0 = bart_thr_fnr0_test_glm, fnr05 = bart_thr_fnr05_test_glm, tss = bart_thr_tss_test_glm),
      bart_test_lm = list(fnr0 = bart_thr_fnr0_test_lm, fnr05 = bart_thr_fnr05_test_lm, tss = bart_thr_tss_test_lm)
    ),
    best_match = best_match,
    overfitting_flags = overfitting_flags,
    moran_summary = if (do_spatial_analysis) {
      data.frame(
        Model = c("bart.step (Train)", best_match, "bart.step (Train vs Test)"),
        k = c(optimal_k_bart$optimal_k, optimal_k_bart_test_lm$optimal_k, optimal_k_bart_test_lm$optimal_k),
        Moran_I = c(moran_bart_train$estimate[1], moran_bart_test_lm$estimate[1], cross_moran_result$estimate[1]),
        P_Value = c(moran_bart_train$p.value, moran_bart_test_lm$p.value, cross_moran_result$p.value)
      )
    } else {
      NULL
    },
    variable_stats = stats_summary,
    best_test_pred = best_test_pred,
    best_test_metrics = list(
      fnr0 = best_test_metrics_fnr0,
      fnr05 = best_test_metrics_fnr05,
      tss = best_test_metrics_tss
    ),
    best_test_thresholds = list(
      fnr0 = best_test_thr_fnr0,
      fnr05 = best_test_thr_fnr05,
      tss = best_test_thr_tss
    )
  )
  
  closeAllConnections()
  while (dev.cur() > 1) dev.off()
  
  return(results)
}
# Example Run
results <- my_bart_step(
  x.data = train_data_normalized, 
  y.data = train_data_normalized$nf, 
  test_data = test_data_normalized,
  selected_vars = c("wc2.1_10m_bio_14", "wc2.1_10m_bio_15", "wc2.1_10m_bio_3","y"  ),#,"x"
  tree.step = 10,
  iter.plot = 5,
  iter.step = 5,
  full = TRUE, 
  quiet = FALSE,
  num_trees = 200, 
  num_burn_in = 250, 
  num_iterations_after_burn_in = 1000, 
  target_fnr = 0.04,
  glm_family = "quasibinomial",
  auc_weight = 0.5,
  ks_weight = 0.5
)
