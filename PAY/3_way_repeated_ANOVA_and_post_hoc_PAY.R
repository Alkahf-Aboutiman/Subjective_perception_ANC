# =========================================================
# Three-way repeated-measures ANOVA on perceived annoyance
# =========================================================

# -------------------------
# Required libraries
# -------------------------
library(conflicted)
library(tidyverse)
library(emmeans)

# Resolve common conflicts explicitly
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

# -------------------------
# Data assumptions
# -------------------------
# data must contain:
# - Pay_Score (numeric, perceived annoyance)
# - Participants (factor)
# - Control_Condition (factor)
# - Noise_Type (factor)
# - Noise_Level (factor)

# -------------------------
# Ensure correct data types
# -------------------------
data <- data %>%
  mutate(
    Participants = as.factor(Participants),
    Control_Condition = as.factor(Control_Condition),
    Noise_Type = as.factor(Noise_Type),
    Noise_Level = as.factor(Noise_Level)
  )

# =========================
# 1. Three-way RM ANOVA
# =========================

anova_pay <- aov(
  Pay_Score ~ Control_Condition * Noise_Type * Noise_Level +
    Error(Participants / (Control_Condition * Noise_Type * Noise_Level)),
  data = data
)

# Display ANOVA results
summary(anova_pay)

# =========================
# 2. Extract F and p-values
# =========================

anova_summary <- summary(anova_pay)

anova_results <- data.frame(
  Effect = character(),
  F_value = numeric(),
  p_value = numeric()
)

# Loop over all ANOVA strata
for (i in seq_along(anova_summary)) {
  
  tab <- anova_summary[[i]][[1]]
  
  if (all(c("F value", "Pr(>F)") %in% colnames(tab))) {
    
    tmp <- data.frame(
      Effect = rownames(tab),
      F_value = tab[, "F value"],
      p_value = tab[, "Pr(>F)"]
    )
    
    anova_results <- rbind(anova_results, tmp)
  }
}

# Remove residual rows
anova_results <- dplyr::filter(anova_results, Effect != "Residuals")

# Apply FDR correction
anova_results$p_value_fdr <- p.adjust(anova_results$p_value, method = "fdr")

print(anova_results)

# =========================
# 3. Tukey post-hoc tests
# =========================

# --- Main effect: Noise Type ---
emm_noise_type <- emmeans(anova_pay, ~ Noise_Type)
tukey_noise_type <- pairs(emm_noise_type, adjust = "tukey")
noise_type_results <- as.data.frame(tukey_noise_type)
print(noise_type_results)

# --- Main effect: Control Condition ---
emm_control <- emmeans(anova_pay, ~ Control_Condition)
tukey_control <- pairs(emm_control, adjust = "tukey")
control_results <- as.data.frame(tukey_control)
print(control_results)

# --- Main effect: Noise Level ---
emm_noise_level <- emmeans(anova_pay, ~ Noise_Level)
tukey_noise_level <- pairs(emm_noise_level, adjust = "tukey")
noise_level_results <- as.data.frame(tukey_noise_level)
print(noise_level_results)

# =========================
# 4. Tukey post-hoc tests for interactions
# =========================

# --- Interaction: Noise Type × Control Condition ---
emm_type_control <- emmeans(anova_pay, ~ Noise_Type * Control_Condition)
tukey_type_control <- pairs(emm_type_control, adjust = "tukey")
type_control_results <- as.data.frame(tukey_type_control)
print(type_control_results)

# --- Interaction: Noise Type × Noise Level ---
emm_type_level <- emmeans(anova_pay, ~ Noise_Type * Noise_Level)
tukey_type_level <- pairs(emm_type_level, adjust = "tukey")
type_level_results <- as.data.frame(tukey_type_level)
print(type_level_results)

# --- Interaction: Control Condition × Noise Level ---
emm_control_level <- emmeans(anova_pay, ~ Control_Condition * Noise_Level)
tukey_control_level <- pairs(emm_control_level, adjust = "tukey")
control_level_results <- as.data.frame(tukey_control_level)
print(control_level_results)

# =========================
# 5. Three-way interaction post-hoc
# =========================

emm_3way <- emmeans(
  anova_pay,
  ~ Control_Condition * Noise_Type * Noise_Level
)

tukey_3way <- pairs(emm_3way, adjust = "tukey")
three_way_results <- as.data.frame(tukey_3way)
print(three_way_results)

# =========================
# 6. Add effect size proxy
# =========================
# estimate corresponds to the mean difference between conditions

three_way_results <- three_way_results %>%
  mutate(abs_diff_mean = abs(estimate))

print(three_way_results)

