required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "car", "glmnet",
  "randomForest", "caret", "xgboost", "e1071", "kernlab"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste(
      "Please install these packages first:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

if (!dir.exists("Model Evaluation")) {
  dir.create("Model Evaluation")
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(glmnet)
library(randomForest)
library(caret)
library(xgboost)
library(e1071)
library(kernlab)

# ---------------------------
# User input
# ---------------------------
data_path <- "ML_mRNA-LNP_FI.csv"

if (!file.exists(data_path)) {
  stop(paste("Data file not found:", data_path))
}

data <- read.csv(data_path, stringsAsFactors = FALSE)

rename_if_present <- function(df, old_name, new_name) {
  if (old_name %in% names(df)) {
    names(df)[names(df) == old_name] <- new_name
  }
  df
}

data <- rename_if_present(data, "Ionizable..1.4.", "Ionizable")
data <- rename_if_present(data, "Phospho..1..2.", "Phospho")
data <- rename_if_present(data, "Helper..1.3.", "Helper")
data <- rename_if_present(data, "PEG..1.2.", "PEG")
data <- rename_if_present(data, "X..amines.in.ionizable.lipid", "Amines_in_Lipid")
data <- rename_if_present(data, "Total.ionizable.amines..nM.", "Total_Amines_nM")
data <- rename_if_present(data, "Size..d.nm.", "Size")
data <- rename_if_present(data, "Zeta..mV.", "Zeta")
data <- rename_if_present(data, "EE.", "EE")

compound_col <- if ("Coumpound.ID" %in% names(data)) {
  "Coumpound.ID"
} else if ("Compound.ID" %in% names(data)) {
  "Compound.ID"
} else {
  stop("Neither 'Coumpound.ID' nor 'Compound.ID' was found in the dataset.")
}

required_columns <- c(
  compound_col,
  "Tumor.expression",
  "Tumor",
  "Ionizable", "Phospho", "Helper", "PEG",
  "Amines_in_Lipid", "Size", "Zeta", "PDI", "EE"
)

missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

data <- data %>%
  group_by(.data[[compound_col]]) %>%
  mutate(rep = row_number()) %>%
  ungroup()

data$Ionizable <- as.factor(data$Ionizable)
data$Phospho  <- as.factor(data$Phospho)
data$Helper   <- as.factor(data$Helper)
data$PEG      <- as.factor(data$PEG)

data$logged_tumor_mrna_expression <- log(data$Tumor.expression)
data$logged_tumor_rhb_biodistribution <- log(data$Tumor)

final_data <- data %>%
  group_by(.data[[compound_col]]) %>%
  summarise(
    logged_tumor_mrna_expression = mean(logged_tumor_mrna_expression, na.rm = TRUE),
    logged_tumor_rhb_biodistribution = mean(logged_tumor_rhb_biodistribution, na.rm = TRUE),
    .groups = "drop"
  )

data_choose <- data[data$rep == 1, ]
data_choose1 <- data_choose %>%
  select(-c(logged_tumor_mrna_expression, logged_tumor_rhb_biodistribution, rep))

data_model <- merge(final_data, data_choose1, by = compound_col)

data_model$Amines_in_Lipid_scale <- as.numeric(scale(data_model$Amines_in_Lipid, scale = TRUE))
data_model$Size_scale <- as.numeric(scale(data_model$Size, scale = TRUE))
data_model$Zeta_scale <- as.numeric(scale(data_model$Zeta, scale = TRUE))
data_model$PDI_scale <- as.numeric(scale(data_model$PDI, scale = TRUE))
data_model$EE_scale <- as.numeric(scale(data_model$EE, scale = TRUE))

empty_metric_table <- function() {
  data.frame(
    Model = c("LM", "RF", "LASSO", "ENET", "XGB", "SVM"),
    RMSE = NA_real_,
    MAE = NA_real_
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# =========================================================
# Physicochemical model evaluation
# =========================================================
calculate_metrics_phys <- function(outcome_name, train_data, test_data) {
  options(warn = -1)
  suppressWarnings({
    formula_str <- paste(
      outcome_name,
      "~ Amines_in_Lipid_scale + Size_scale + Zeta_scale + PDI_scale + EE_scale"
    )
    
    y_train <- train_data[[outcome_name]]
    
    if (length(unique(y_train[!is.na(y_train)])) < 2) {
      return(empty_metric_table())
    }
    
    lm_model <- tryCatch(
      lm(as.formula(formula_str), data = train_data),
      error = function(e) NULL
    )
    
    X_train <- tryCatch(
      model.matrix(as.formula(formula_str), data = train_data)[, -1, drop = FALSE],
      error = function(e) NULL
    )
    X_test <- tryCatch(
      model.matrix(as.formula(formula_str), data = test_data)[, -1, drop = FALSE],
      error = function(e) NULL
    )
    
    if (is.null(X_train) || is.null(X_test) || nrow(X_train) == 0 || ncol(X_train) == 0) {
      return(empty_metric_table())
    }
    
    # LM
    lm_pred <- tryCatch(
      as.numeric(predict(lm_model, test_data)),
      error = function(e) rep(NA, nrow(test_data))
    )
    
    # LASSO
    lasso_pred <- tryCatch({
      set.seed(123)
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian")
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = cv_model$lambda.min)
      as.numeric(predict(lasso_model, X_test))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # Elastic Net
    enet_pred <- tryCatch({
      alphas <- seq(0, 1, by = 0.1)
      cv_results <- list()
      cv_errors <- numeric(length(alphas))
      for (i in seq_along(alphas)) {
        cv_fit <- cv.glmnet(X_train, y_train, alpha = alphas[i], family = "gaussian")
        cv_results[[i]] <- cv_fit
        cv_errors[i] <- min(cv_fit$cvm)
      }
      best_alpha <- alphas[which.min(cv_errors)]
      best_model <- cv_results[[which.min(cv_errors)]]
      enet_model <- glmnet(X_train, y_train, alpha = best_alpha, lambda = best_model$lambda.min)
      as.numeric(predict(enet_model, X_test))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # RF
    rf_pred <- tryCatch({
      train_control <- trainControl(method = "cv", number = 5)
      grid <- expand.grid(mtry = c(3, 4, 5))
      set.seed(111)
      forest_gridsearch <- train(
        as.formula(formula_str),
        data = train_data,
        method = "rf",
        trControl = train_control,
        tuneGrid = grid,
        ntree = 500
      )
      best_mtry <- forest_gridsearch$bestTune$mtry
      set.seed(111)
      forest_model <- randomForest(
        as.formula(formula_str),
        data = train_data,
        mtry = best_mtry,
        ntree = 500,
        importance = TRUE
      )
      as.numeric(predict(forest_model, test_data))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # XGB
    xgb_pred <- tryCatch({
      xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
      xgb_grid <- expand.grid(
        nrounds = c(100, 200),
        max_depth = c(3, 6),
        eta = c(0.1, 0.3),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
      )
      set.seed(123)
      xgb_gridsearch <- train(
        x = X_train, y = y_train,
        method = "xgbTree",
        trControl = trainControl(method = "cv", number = 3),
        tuneGrid = xgb_grid,
        verbosity = 0
      )
      best_params <- xgb_gridsearch$bestTune
      params <- list(
        objective = "reg:squarederror",
        max_depth = best_params$max_depth,
        eta = best_params$eta,
        gamma = best_params$gamma,
        colsample_bytree = best_params$colsample_bytree,
        min_child_weight = best_params$min_child_weight,
        subsample = best_params$subsample
      )
      set.seed(123)
      xgb_model <- xgb.train(
        params = params,
        data = xgb_train,
        nrounds = best_params$nrounds,
        verbose = 0
      )
      as.numeric(predict(xgb_model, as.matrix(X_test)))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # SVM
    svm_pred <- tryCatch({
      set.seed(123)
      svm_model <- train(
        as.formula(formula_str),
        data = train_data,
        method = "svmLinear",
        trControl = trainControl(method = "cv", number = 3),
        tuneGrid = expand.grid(C = c(0.1, 1, 10))
      )
      as.numeric(predict(svm_model, test_data))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    actual <- test_data[[outcome_name]]
    
    RMSE_vec <- c(
      sqrt(mean((actual - lm_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - rf_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - lasso_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - enet_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - xgb_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - svm_pred)^2, na.rm = TRUE))
    )
    
    MAE_vec <- c(
      mean(abs(actual - lm_pred), na.rm = TRUE),
      mean(abs(actual - rf_pred), na.rm = TRUE),
      mean(abs(actual - lasso_pred), na.rm = TRUE),
      mean(abs(actual - enet_pred), na.rm = TRUE),
      mean(abs(actual - xgb_pred), na.rm = TRUE),
      mean(abs(actual - svm_pred), na.rm = TRUE)
    )
    
    data.frame(
      Model = c("LM", "RF", "LASSO", "ENET", "XGB", "SVM"),
      RMSE = round(RMSE_vec, 3),
      MAE = round(MAE_vec, 3)
    )
  })
}

perform_cv_phys <- function(df, outcomes, n_folds = 5) {
  set.seed(123)
  compound_ids <- unique(df[[compound_col]])
  n_folds <- min(n_folds, length(compound_ids))
  compounds_folds <- data.frame(
    compound_id = compound_ids,
    fold = sample(rep(1:n_folds, length.out = length(compound_ids)))
  )
  
  all_cv_results <- list()
  
  for (outcome_name in outcomes) {
    cv_metrics <- list()
    
    for (fold in 1:n_folds) {
      test_compounds <- compounds_folds$compound_id[compounds_folds$fold == fold]
      test_data <- df[df[[compound_col]] %in% test_compounds, ]
      train_data <- df[!df[[compound_col]] %in% test_compounds, ]
      cv_metrics[[fold]] <- calculate_metrics_phys(outcome_name, train_data, test_data)
    }
    
    avg_metrics <- do.call(rbind, cv_metrics) %>%
      group_by(Model) %>%
      summarise(
        RMSE = safe_mean(RMSE),
        MAE = safe_mean(MAE),
        .groups = "drop"
      ) %>%
      mutate(Outcome = outcome_name, .before = 1)
    
    all_cv_results[[outcome_name]] <- avg_metrics
  }
  
  bind_rows(all_cv_results)
}

# =========================================================
# Composition model evaluation
# =========================================================
calculate_metrics_comp <- function(outcome_name, train_data, test_data) {
  options(warn = -1)
  suppressWarnings({
    formula_str <- paste(outcome_name, "~ Ionizable + Phospho + Helper + PEG")
    
    y_train <- train_data[[outcome_name]]
    
    if (length(unique(y_train[!is.na(y_train)])) < 2) {
      return(empty_metric_table())
    }
    
    lm_model <- tryCatch(
      lm(as.formula(formula_str), data = train_data),
      error = function(e) NULL
    )
    
    X_train <- tryCatch(
      model.matrix(as.formula(formula_str), data = train_data)[, -1, drop = FALSE],
      error = function(e) NULL
    )
    X_test <- tryCatch(
      model.matrix(as.formula(formula_str), data = test_data)[, -1, drop = FALSE],
      error = function(e) NULL
    )
    
    if (is.null(X_train) || is.null(X_test) || nrow(X_train) == 0 || ncol(X_train) == 0) {
      return(empty_metric_table())
    }
    
    # LM
    lm_pred <- tryCatch(
      as.numeric(predict(lm_model, test_data)),
      error = function(e) rep(NA, nrow(test_data))
    )
    
    # LASSO
    lasso_pred <- tryCatch({
      set.seed(123)
      cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian")
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = cv_lasso$lambda.min)
      as.numeric(predict(lasso_model, X_test))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # ENET
    enet_pred <- tryCatch({
      alphas <- seq(0, 1, by = 0.1)
      cv_errors <- numeric(length(alphas))
      cv_results <- vector("list", length(alphas))
      for (i in seq_along(alphas)) {
        cv_fit <- cv.glmnet(X_train, y_train, alpha = alphas[i], family = "gaussian")
        cv_results[[i]] <- cv_fit
        cv_errors[i] <- min(cv_fit$cvm)
      }
      best_alpha <- alphas[which.min(cv_errors)]
      best_model <- cv_results[[which.min(cv_errors)]]
      enet_model <- glmnet(X_train, y_train, alpha = best_alpha, lambda = best_model$lambda.min)
      as.numeric(predict(enet_model, X_test))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # RF
    rf_pred <- tryCatch({
      set.seed(111)
      rf_grid <- expand.grid(mtry = c(3, 4, 5))
      rf_cv <- trainControl(method = "cv", number = 5)
      rf_tune <- train(
        as.formula(formula_str),
        data = train_data,
        method = "rf",
        trControl = rf_cv,
        tuneGrid = rf_grid,
        ntree = 500
      )
      best_mtry <- rf_tune$bestTune$mtry
      forest_model <- randomForest(
        as.formula(formula_str),
        data = train_data,
        mtry = best_mtry,
        ntree = 500,
        importance = TRUE
      )
      as.numeric(predict(forest_model, test_data))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # XGB
    xgb_pred <- tryCatch({
      xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
      xgb_grid <- expand.grid(
        nrounds = c(100, 200),
        max_depth = c(3, 6),
        eta = c(0.1, 0.3),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
      )
      set.seed(123)
      xgb_tune <- train(
        x = X_train, y = y_train,
        method = "xgbTree",
        trControl = trainControl(method = "cv", number = 3),
        tuneGrid = xgb_grid,
        verbosity = 0
      )
      best_params <- xgb_tune$bestTune
      params <- list(
        objective = "reg:squarederror",
        max_depth = best_params$max_depth,
        eta = best_params$eta,
        gamma = best_params$gamma,
        colsample_bytree = best_params$colsample_bytree,
        min_child_weight = best_params$min_child_weight,
        subsample = best_params$subsample
      )
      set.seed(123)
      xgb_model <- xgb.train(
        params = params,
        data = xgb_train,
        nrounds = best_params$nrounds,
        verbose = 0
      )
      as.numeric(predict(xgb_model, as.matrix(X_test)))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    # SVM
    svm_pred <- tryCatch({
      set.seed(123)
      svm_model <- train(
        x = X_train,
        y = y_train,
        method = "svmLinear",
        trControl = trainControl(method = "cv", number = 3),
        tuneGrid = expand.grid(C = c(0.1, 1, 10))
      )
      as.numeric(predict(svm_model, X_test))
    }, error = function(e) rep(NA, nrow(test_data)))
    
    actual <- test_data[[outcome_name]]
    
    RMSE_vec <- c(
      sqrt(mean((actual - lm_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - rf_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - lasso_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - enet_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - xgb_pred)^2, na.rm = TRUE)),
      sqrt(mean((actual - svm_pred)^2, na.rm = TRUE))
    )
    
    MAE_vec <- c(
      mean(abs(actual - lm_pred), na.rm = TRUE),
      mean(abs(actual - rf_pred), na.rm = TRUE),
      mean(abs(actual - lasso_pred), na.rm = TRUE),
      mean(abs(actual - enet_pred), na.rm = TRUE),
      mean(abs(actual - xgb_pred), na.rm = TRUE),
      mean(abs(actual - svm_pred), na.rm = TRUE)
    )
    
    data.frame(
      Model = c("LM", "RF", "LASSO", "ENET", "XGB", "SVM"),
      RMSE = round(RMSE_vec, 3),
      MAE = round(MAE_vec, 3)
    )
  })
}

perform_cv_comp <- function(df, outcomes, n_folds = 5) {
  set.seed(123)
  compound_ids <- unique(df[[compound_col]])
  n_folds <- min(n_folds, length(compound_ids))
  compounds_folds <- data.frame(
    compound_id = compound_ids,
    fold = sample(rep(1:n_folds, length.out = length(compound_ids)))
  )
  
  all_cv_results <- list()
  
  for (outcome_name in outcomes) {
    cv_metrics <- list()
    
    for (fold in 1:n_folds) {
      test_compounds <- compounds_folds$compound_id[compounds_folds$fold == fold]
      test_data <- df[df[[compound_col]] %in% test_compounds, ]
      train_data <- df[!df[[compound_col]] %in% test_compounds, ]
      cv_metrics[[fold]] <- calculate_metrics_comp(outcome_name, train_data, test_data)
    }
    
    avg_metrics <- do.call(rbind, cv_metrics) %>%
      group_by(Model) %>%
      summarise(
        RMSE = safe_mean(RMSE),
        MAE = safe_mean(MAE),
        .groups = "drop"
      ) %>%
      mutate(Outcome = outcome_name, .before = 1)
    
    all_cv_results[[outcome_name]] <- avg_metrics
  }
  
  bind_rows(all_cv_results)
}

# =========================================================
# Run evaluations
# =========================================================
phys_eval <- perform_cv_phys(
  data_model,
  c("logged_tumor_mrna_expression", "logged_tumor_rhb_biodistribution"),
  n_folds = 5
)

phys_eval$Outcome <- dplyr::case_when(
  phys_eval$Outcome == "logged_tumor_mrna_expression" ~ "logged tumor in mRNA expression",
  phys_eval$Outcome == "logged_tumor_rhb_biodistribution" ~ "logged tumor in RhB biodistribution",
  TRUE ~ phys_eval$Outcome
)

write.csv(
  phys_eval,
  "Model Evaluation/physicochemical_model_evaluation_metrics.csv",
  row.names = FALSE
)

comp_eval <- perform_cv_comp(
  data_model,
  c("logged_tumor_mrna_expression", "logged_tumor_rhb_biodistribution"),
  n_folds = 5
)

comp_eval$Outcome <- dplyr::case_when(
  comp_eval$Outcome == "logged_tumor_mrna_expression" ~ "logged tumor in mRNA expression",
  comp_eval$Outcome == "logged_tumor_rhb_biodistribution" ~ "logged tumor in RhB biodistribution",
  TRUE ~ comp_eval$Outcome
)

write.csv(
  comp_eval,
  "Model Evaluation/composition_model_evaluation_metrics.csv",
  row.names = FALSE
)

message("Done.")
message("Saved files:")
message("Model Evaluation/physicochemical_model_evaluation_metrics.csv")
message("Model Evaluation/composition_model_evaluation_metrics.csv")