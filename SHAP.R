library(Matrix)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(MuMIn)
library(kableExtra)
library(knitr)
library(randomForest)
library(glmnet)
library(car)
library(xgboost)
library(iml)
library(scales)
library(patchwork)

data = read.csv("ML_mRNA-LNP_FI - Sample.csv")

data <- data %>%
  group_by(Coumpound.ID) %>%
  mutate(rep = row_number()) %>%
  ungroup()

names(data)[which(names(data)=="Ionizable..1.4.")] = "Ionizable"
names(data)[which(names(data)=="Phospho..1..2.")] = "Phospho"
names(data)[which(names(data)=="Helper..1.3.")] = "Helper"
names(data)[which(names(data)=="PEG..1.2.")] = "PEG"

data$Ionizable = as.factor(data$Ionizable)
data$Phospho = as.factor(data$Phospho)
data$Helper = as.factor(data$Helper)
data$PEG = as.factor(data$PEG)

names(data)[which(names(data)=="X..amines.in.ionizable.lipid")] = "Amines_in_Lipid"
names(data)[which(names(data)=="Total.ionizable.amines..nM.")] = "Total_Amines_nM"
names(data)[which(names(data)=="Size..d.nm.")] = "Size"
names(data)[which(names(data)=="Zeta..mV.")] = "Zeta"
names(data)[which(names(data)=="PDI")] = "PDI"
names(data)[which(names(data)=="EE.")] = "EE"

y0 = log(data$Tumor.expression)
y1 = log(data$Tumor)

data = cbind(data, y0, y1)

final_data <- data %>%
  group_by(Coumpound.ID) %>%
  summarize(across(y0:y1, mean, na.rm = TRUE))

data_choose = data[data$rep == 1, ]
data_choose1 = data_choose %>% select(-c(y0, y1, rep))
data = merge(final_data, data_choose1, by="Coumpound.ID")

data$Amines_in_Lipid_scale = as.numeric(scale(data$Amines_in_Lipid))
data$Size_scale = as.numeric(scale(data$Size))
data$Zeta_scale = as.numeric(scale(data$Zeta))
data$PDI_scale = as.numeric(scale(data$PDI))
data$EE_scale = as.numeric(scale(data$EE))

if (!dir.exists("SHAP")) {
  dir.create("SHAP")
}

# ============================
# Common plot theme
# ============================
plot_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.4),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

# ============================
# 1. PHYSICOCHEMICAL SHAP
# ============================

set.seed(123)

predictors <- c("Amines_in_Lipid_scale","Size_scale","Zeta_scale","PDI_scale","EE_scale")
B <- 50
outcome_name <- "y1"

rf_model_ref <- randomForest(
  as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + "))),
  data = data,
  ntree = 500
)

predictor_rf_ref <- Predictor$new(
  rf_model_ref,
  data = data[, predictors],
  y = data[[outcome_name]]
)

shapley_list_ref <- lapply(1:nrow(data), function(i) {
  Shapley$new(predictor_rf_ref, x.interest = data[i, predictors, drop = FALSE])$results
})

shapley_all_ref <- do.call(rbind, shapley_list_ref)

shapley_ref_summary <- shapley_all_ref %>%
  group_by(feature) %>%
  summarise(mean_phi_ref = mean(phi, na.rm = TRUE), .groups = "drop")

boot_results <- list()

for (b in 1:B) {
  set.seed(1000 + b)
  boot_idx <- sample(1:nrow(data), replace = TRUE)
  boot_data <- data[boot_idx, ]
  
  rf_model <- randomForest(
    as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + "))),
    data = boot_data,
    ntree = 500
  )
  
  predictor_rf <- Predictor$new(
    rf_model,
    data = boot_data[, predictors],
    y = boot_data[[outcome_name]]
  )
  
  shapley_list <- lapply(1:nrow(boot_data), function(i) {
    Shapley$new(predictor_rf, x.interest = boot_data[i, predictors, drop = FALSE])$results
  })
  
  shapley_all <- do.call(rbind, shapley_list)
  
  boot_results[[b]] <- shapley_all %>%
    group_by(feature) %>%
    summarise(mean_phi = mean(phi, na.rm = TRUE), .groups = "drop")
}

boot_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), boot_results)
colnames(boot_df)[-1] <- paste0("boot", 1:B)

boot_summary_phys <- boot_df %>%
  pivot_longer(-feature, names_to = "boot", values_to = "phi") %>%
  group_by(feature) %>%
  summarise(
    se = sd(phi, na.rm = TRUE) / sqrt(B),
    mean_phi_ref = shapley_ref_summary$mean_phi_ref[match(feature, shapley_ref_summary$feature)],
    lower = mean_phi_ref - 1.96 * se,
    upper = mean_phi_ref + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_phi_ref > 0, "Positive", "Negative"))

p1 <- ggplot(boot_summary_phys,
             aes(x = feature, y = mean_phi_ref,
                 color = direction, shape = direction)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.6, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.22, linewidth = 0.8) +
  geom_point(size = 4.5) +
  coord_flip() +
  scale_color_manual(values = c("Negative" = "#377EB8", "Positive" = "#E41A1C")) +
  scale_shape_manual(values = c("Negative" = 16, "Positive" = 17)) +
  labs(title = "Shapley Values – Physicochemical",
       x = "Feature",
       y = "Shapley value") +
  plot_theme

# ============================
# 2. COMPOSITION SHAP
# ============================

predictors <- c("Ionizable","Phospho","Helper","PEG")

boot_results <- list()

for (b in 1:B) {
  set.seed(123 + b)
  boot_idx <- sample(1:nrow(data), replace = TRUE)
  boot_data <- data[boot_idx, ]
  
  rf_model <- randomForest(
    as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + "))),
    data = boot_data,
    ntree = 500
  )
  
  shap_data <- boot_data[, predictors]
  
  predictor_rf <- Predictor$new(rf_model, data = shap_data, y = boot_data[[outcome_name]])
  
  shapley_list <- lapply(1:nrow(shap_data), function(i) {
    Shapley$new(predictor_rf, x.interest = shap_data[i, , drop = FALSE])$results
  })
  
  shapley_all <- do.call(rbind, shapley_list)
  
  shapley_signed <- shapley_all %>%
    group_by(feature) %>%
    summarise(mean_phi = mean(phi, na.rm = TRUE), .groups = "drop")
  
  boot_results[[b]] <- shapley_signed
}

boot_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), boot_results)
colnames(boot_df)[-1] <- paste0("boot", 1:B)

boot_summary_comp <- boot_df %>%
  pivot_longer(-feature, names_to = "boot", values_to = "phi") %>%
  group_by(feature) %>%
  summarise(
    mean_phi = mean(phi, na.rm = TRUE),
    se = sd(phi, na.rm = TRUE) / sqrt(B),
    lower = mean_phi - 1.96 * se,
    upper = mean_phi + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_phi > 0, "Positive", "Negative"))

p2 <- ggplot(boot_summary_comp,
             aes(x = feature, y = mean_phi,
                 color = direction, shape = direction)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.6, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.22, linewidth = 0.8) +
  geom_point(size = 4.5) +
  coord_flip() +
  scale_color_manual(values = c("Negative" = "#377EB8", "Positive" = "#E41A1C")) +
  scale_shape_manual(values = c("Negative" = 16, "Positive" = 17)) +
  labs(title = "Shapley Values – Composition",
       x = "Feature",
       y = "Shapley value") +
  plot_theme

# ============================
# COMBINED FIGURE
# ============================

combined_plot <- p2 + p1 + plot_layout(ncol = 2)

ggsave("SHAP/shap_combined.png",
       combined_plot,
       width = 16, height = 6.5,
       dpi = 600, bg = "white")

print(combined_plot)