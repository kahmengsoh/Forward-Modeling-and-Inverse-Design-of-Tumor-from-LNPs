# Inverse design using Random Forest + MCMC
# Clean R script for GitHub / paper reference

# ---------------------------
# Required packages
# ---------------------------
required_packages <- c("dplyr", "tidyr", "caret", "randomForest", "ggplot2")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste(
      "Please install the following packages before running this script:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(dplyr)
library(tidyr)
library(caret)
library(randomForest)
library(ggplot2)

# ---------------------------
# User input
# ---------------------------
data_path <- "ML_mRNA-LNP_FI.csv"
output_dir <- "MCMC"

if (!file.exists(data_path)) {
  stop(paste("Data file not found:", data_path))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

set.seed(123)

# ---------------------------
# Load data
# ---------------------------
data <- read.csv(data_path, stringsAsFactors = FALSE)

# ---------------------------
# Rename variables if needed
# ---------------------------
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

for (v in c("Ionizable", "Phospho", "Helper", "PEG")) {
  if (v %in% names(data)) {
    data[[v]] <- as.factor(data[[v]])
  }
}

# ---------------------------
# Check required columns
# ---------------------------
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
  "Amines_in_Lipid",
  "Size",
  "Zeta",
  "PDI",
  "EE"
)

missing_columns <- setdiff(required_columns, names(data))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required columns are missing from the dataset:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ---------------------------
# Add replicate index
# ---------------------------
data <- data %>%
  group_by(.data[[compound_col]]) %>%
  mutate(rep = row_number()) %>%
  ungroup()

# ---------------------------
# Create outcomes
# ---------------------------
data$tumor_mrna_expression <- log(data$Tumor.expression)
data$tumor_rhb_biodistribution <- log(data$Tumor)

# ---------------------------
# Aggregate outcomes by compound
# ---------------------------
final_data <- data %>%
  group_by(.data[[compound_col]]) %>%
  summarise(
    tumor_mrna_expression = mean(tumor_mrna_expression, na.rm = TRUE),
    tumor_rhb_biodistribution = mean(tumor_rhb_biodistribution, na.rm = TRUE),
    .groups = "drop"
  )

# Keep one row of predictors per compound
data_choose <- data[data$rep == 1, ]
data_choose1 <- data_choose %>%
  select(-c(tumor_mrna_expression, tumor_rhb_biodistribution, rep))

data_model <- merge(final_data, data_choose1, by = compound_col)

# ---------------------------
# Helper functions
# ---------------------------
log_posterior <- function(X, model, target_y, observed_data) {
  x_df <- data.frame(
    Amines_in_Lipid = X[1],
    Size = X[2],
    Zeta = X[3],
    PDI = X[4],
    EE = X[5]
  )
  
  X_bounds <- list(
    Amines_in_Lipid = range(observed_data$Amines_in_Lipid, na.rm = TRUE),
    Size = range(observed_data$Size, na.rm = TRUE),
    Zeta = range(observed_data$Zeta, na.rm = TRUE),
    PDI = range(observed_data$PDI, na.rm = TRUE),
    EE = range(observed_data$EE, na.rm = TRUE)
  )
  
  within_bounds <- all(
    X[1] >= X_bounds$Amines_in_Lipid[1] & X[1] <= X_bounds$Amines_in_Lipid[2],
    X[2] >= X_bounds$Size[1] & X[2] <= X_bounds$Size[2],
    X[3] >= X_bounds$Zeta[1] & X[3] <= X_bounds$Zeta[2],
    X[4] >= X_bounds$PDI[1] & X[4] <= X_bounds$PDI[2],
    X[5] >= X_bounds$EE[1] & X[5] <= X_bounds$EE[2]
  )
  
  if (!within_bounds) return(-Inf)
  
  pred_y <- predict(model, newdata = x_df)
  sigma <- 0.15
  dnorm(pred_y, mean = target_y, sd = sigma, log = TRUE)
}

run_metropolis <- function(log_post_func, init_values, n_iter, proposal_sd, model, target_y, observed_data) {
  n_params <- length(init_values)
  chain <- matrix(0, nrow = n_iter, ncol = n_params)
  chain[1, ] <- init_values
  
  for (i in 2:n_iter) {
    proposal <- chain[i - 1, ] + rnorm(n_params, 0, proposal_sd)
    
    log_alpha <- log_post_func(proposal, model, target_y, observed_data) -
      log_post_func(chain[i - 1, ], model, target_y, observed_data)
    
    if (log(runif(1)) < log_alpha) {
      chain[i, ] <- proposal
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  
  chain
}

# ---------------------------
# Random Forest + inverse design
# ---------------------------
train_control <- trainControl(method = "cv", number = 5)
rf_results_list <- list()
n_boot <- 50

for (outcome_name in c("tumor_mrna_expression", "tumor_rhb_biodistribution")) {
  formula_str <- paste(outcome_name, "~ Amines_in_Lipid + Size + Zeta + PDI + EE")
  grid <- expand.grid(mtry = c(3, 4, 5))
  
  set.seed(123)
  rf_gridsearch <- train(
    as.formula(formula_str),
    data = data_model,
    method = "rf",
    trControl = train_control,
    tuneGrid = grid,
    ntree = 500
  )
  
  best_mtry <- rf_gridsearch$bestTune$mtry
  
  set.seed(123)
  model <- randomForest(
    as.formula(formula_str),
    data = data_model,
    mtry = best_mtry,
    ntree = 500,
    importance = TRUE
  )
  
  multiples <- c(0.7, 0.8, 0.9, 1.0)
  max_val <- max(data_model[[outcome_name]], na.rm = TRUE)
  target_Y_values <- multiples * max_val
  target_labels <- paste0(multiples, "xMax")
  
  single_run_list <- list()
  
  for (i in seq_along(target_Y_values)) {
    target_Y <- target_Y_values[i]
    target_label <- target_labels[i]
    
    initial_X <- c(
      median(data_model$Amines_in_Lipid, na.rm = TRUE),
      median(data_model$Size, na.rm = TRUE),
      median(data_model$Zeta, na.rm = TRUE),
      median(data_model$PDI, na.rm = TRUE),
      median(data_model$EE, na.rm = TRUE)
    )
    
    proposal_sd <- c(
      diff(range(data_model$Amines_in_Lipid, na.rm = TRUE)) * 0.05,
      diff(range(data_model$Size, na.rm = TRUE)) * 0.05,
      diff(range(data_model$Zeta, na.rm = TRUE)) * 0.05,
      diff(range(data_model$PDI, na.rm = TRUE)) * 0.05,
      diff(range(data_model$EE, na.rm = TRUE)) * 0.05
    )
    
    set.seed(123)
    chain <- run_metropolis(
      log_post_func = log_posterior,
      init_values = initial_X,
      n_iter = 2000,
      proposal_sd = proposal_sd,
      model = model,
      target_y = target_Y,
      observed_data = data_model
    )
    
    final_chain <- chain[seq(501, 2000, by = 10), ]
    df_chain <- as.data.frame(final_chain)
    colnames(df_chain) <- c("Amines_in_Lipid", "Size", "Zeta", "PDI", "EE")
    
    single_run_list[[target_label]] <- df_chain %>%
      summarise(across(everything(), mean))
  }
  
  single_run_df <- bind_rows(single_run_list, .id = "Target_Multiplier") %>%
    mutate(
      Outcome_Var = outcome_name,
      Target_Y = target_Y_values
    )
  
  # ---------------------------
  # Bootstrap for SE and 95% CI
  # ---------------------------
  boot_list <- list()
  
  for (b in 1:n_boot) {
    set.seed(100 + b)
    boot_data <- data_model %>% sample_frac(replace = TRUE)
    
    model_b <- randomForest(
      as.formula(formula_str),
      data = boot_data,
      mtry = best_mtry,
      ntree = 300
    )
    
    boot_pred <- single_run_df %>%
      rowwise() %>%
      mutate(
        pred_boot = predict(
          model_b,
          newdata = data.frame(
            Amines_in_Lipid = Amines_in_Lipid,
            Size = Size,
            Zeta = Zeta,
            PDI = PDI,
            EE = EE
          )
        )
      ) %>%
      ungroup()
    
    boot_list[[b]] <- boot_pred
  }
  
  boot_all <- bind_rows(boot_list, .id = "boot_id")
  
  se_df <- boot_all %>%
    group_by(Outcome_Var, Target_Multiplier) %>%
    summarise(SE = sd(pred_boot), .groups = "drop")
  
  rf_final <- single_run_df %>%
    left_join(se_df, by = c("Outcome_Var", "Target_Multiplier")) %>%
    mutate(
      across(c(Amines_in_Lipid, Size, Zeta, PDI, EE), ~ .x, .names = "{.col}_Mean")
    ) %>%
    mutate(
      across(c(Amines_in_Lipid, Size, Zeta, PDI, EE), ~ .x - 1.96 * SE, .names = "{.col}_Low"),
      across(c(Amines_in_Lipid, Size, Zeta, PDI, EE), ~ .x + 1.96 * SE, .names = "{.col}_High")
    )
  
  rf_results_list[[outcome_name]] <- rf_final
}

# ---------------------------
# Combine results
# ---------------------------
rf_results_raw <- bind_rows(rf_results_list)

# ---------------------------
# Paper-style tables with CI
# ---------------------------
rf_results_tbl_mrna <- rf_results_raw %>%
  filter(Outcome_Var == "tumor_mrna_expression") %>%
  transmute(
    `Target Logged Tumor in mRNA Expression` = round(Target_Y, 3),
    `Target Multiplier` = Target_Multiplier,
    `Amines in Lipid (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Amines_in_Lipid_Mean, Amines_in_Lipid_Low, Amines_in_Lipid_High),
    `Size (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Size_Mean, Size_Low, Size_High),
    `Zeta (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Zeta_Mean, Zeta_Low, Zeta_High),
    `PDI (95% CI)` = sprintf("%.3f (%.3f-%.3f)", PDI_Mean, PDI_Low, PDI_High),
    `EE (95% CI)` = sprintf("%.3f (%.3f-%.3f)", EE_Mean, EE_Low, EE_High)
  )

rf_results_tbl_rhb <- rf_results_raw %>%
  filter(Outcome_Var == "tumor_rhb_biodistribution") %>%
  transmute(
    `Target Logged Tumor in RhB Biodistribution` = round(Target_Y, 3),
    `Target Multiplier` = Target_Multiplier,
    `Amines in Lipid (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Amines_in_Lipid_Mean, Amines_in_Lipid_Low, Amines_in_Lipid_High),
    `Size (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Size_Mean, Size_Low, Size_High),
    `Zeta (95% CI)` = sprintf("%.3f (%.3f-%.3f)", Zeta_Mean, Zeta_Low, Zeta_High),
    `PDI (95% CI)` = sprintf("%.3f (%.3f-%.3f)", PDI_Mean, PDI_Low, PDI_High),
    `EE (95% CI)` = sprintf("%.3f (%.3f-%.3f)", EE_Mean, EE_Low, EE_High)
  )

# ---------------------------
# Save CSV files
# ---------------------------
write.csv(
  rf_results_tbl_mrna,
  file.path(output_dir, "rf_inverse_design_tumor_mrna_expression.csv"),
  row.names = FALSE
)

write.csv(
  rf_results_tbl_rhb,
  file.path(output_dir, "rf_inverse_design_tumor_rhb_biodistribution.csv"),
  row.names = FALSE
)


# ---------------------------
# Prepare plot data
# ---------------------------
rf_plot_df <- rf_results_raw %>%
  mutate(Multiple = as.numeric(sub("xMax", "", Target_Multiplier))) %>%
  filter(!is.na(Multiple)) %>%
  pivot_longer(
    cols = c(
      Amines_in_Lipid_Mean, Size_Mean, Zeta_Mean, PDI_Mean, EE_Mean,
      Amines_in_Lipid_Low, Size_Low, Zeta_Low, PDI_Low, EE_Low,
      Amines_in_Lipid_High, Size_High, Zeta_High, PDI_High, EE_High
    ),
    names_to = c("Property", ".value"),
    names_pattern = "(.*)_(Mean|Low|High)"
  )

rf_plot_df$Property[rf_plot_df$Property == "Amines_in_Lipid"] <- "Amines in Lipid"

# ---------------------------
# Save plots with CI error bars
# ---------------------------
for (outcome_name in unique(rf_plot_df$Outcome_Var)) {
  df_subset <- filter(rf_plot_df, Outcome_Var == outcome_name)
  
  plot_title <- ifelse(
    outcome_name == "tumor_mrna_expression",
    "RF + Bootstrap CI for tumor in mRNA expression",
    "RF + Bootstrap CI for tumor in RhB biodistribution"
  )
  
  file_name <- ifelse(
    outcome_name == "tumor_mrna_expression",
    "tumor_mrna_expression.png",
    "tumor_rhb_biodistribution.png"
  )
  
  p <- ggplot(df_subset, aes(x = Multiple, y = Mean, color = Property, group = Property)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Low, ymax = High), width = 0.05) +
    facet_wrap(~Property, scales = "free_y") +
    labs(
      title = plot_title,
      x = "Target Y (multiples of Max)",
      y = "Predictor Value"
    ) + theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  ggsave(
    filename = file.path(output_dir, file_name),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

message("Analysis complete.")
message("Saved files in: ", output_dir)