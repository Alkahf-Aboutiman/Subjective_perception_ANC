# =========================================================
# Three-way repeated-measures ANOVA on perceived annoyance
# =========================================================

# -------------------------
# Required libraries
# -------------------------
library(conflicted)
library(tidyverse)
library(emmeans)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)

# -------------------------
# 1. Import CSV data
# -------------------------

# Read PAY scores
pay_data <- read.csv(
  "PAY_results.csv",
  stringsAsFactors = TRUE,
  na.strings = c("", "NA")
)

# Read stimulus metadata
stimuli_data <- read.csv(
  "Stimuli_information.csv",
  stringsAsFactors = TRUE
)

# Replace missing PAY values by neutral value (scale midpoint)
pay_data[is.na(pay_data)] <- 3

# -------------------------
# 2. Reshape PAY data
# -------------------------

# Harmonize stimulus column names
pay_data <- pay_data %>%
  rename_with(~ gsub("Stimulus.", "Stimuli_", .), starts_with("Stimulus"))

# Convert PAY scores to long format
pay_scores <- pay_data %>%
  select(-Participants) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Stimulus",
    values_to = "Pay_Score"
  )

# Prepare stimulus metadata
stimuli_info <- stimuli_data %>%
  rename(Stimulus = Stimuli) %>%
  mutate(Stimulus = as.character(Stimulus))

# Number of stimuli per participant
n_stimuli <- ncol(pay_data) - 1

# Merge PAY scores with stimulus metadata
data <- pay_scores %>%
  left_join(stimuli_info, by = "Stimulus") %>%
  mutate(
    Participants = rep(pay_data$Participants, each = n_stimuli),
    Pay_Score = replace_na(Pay_Score, 3)
  )

# -------------------------
# 3. Ensure correct data types
# -------------------------

data <- data %>%
  mutate(
    Participants = as.factor(Participants),
    Control_Condition = as.factor(Control_Condition),
    Noise_Type = as.factor(Noise_Type),
    Noise_Level = as.factor(Noise_Level)
  )

# =========================
# 4. Three-way RM ANOVA
# =========================

anova_pay <- aov(
  Pay_Score ~ Control_Condition * Noise_Type * Noise_Level +
    Error(Participants / (Control_Condition * Noise_Type * Noise_Level)),
  data = data
)

summary(anova_pay)

# =========================
# 5. Extract F and p-values
# =========================

anova_summary <- summary(anova_pay)

anova_results <- data.frame(
  Effect = character(),
  F_value = numeric(),
  p_value = numeric()
)

for (i in seq_along(anova_summary)) {
  tab <- anova_summary[[i]][[1]]
  
  if (all(c("F value", "Pr(>F)") %in% colnames(tab))) {
    anova_results <- rbind(
      anova_results,
      data.frame(
        Effect = rownames(tab),
        F_value = tab[, "F value"],
        p_value = tab[, "Pr(>F)"]
      )
    )
  }
}

anova_results <- dplyr::filter(anova_results, Effect != "Residuals")
anova_results$p_value_fdr <- p.adjust(anova_results$p_value, method = "fdr")

print(anova_results)

# =========================
# 6. Tukey post-hoc tests
# =========================

# Main effects
pairs(emmeans(anova_pay, ~ Noise_Type), adjust = "tukey")
pairs(emmeans(anova_pay, ~ Control_Condition), adjust = "tukey")
pairs(emmeans(anova_pay, ~ Noise_Level), adjust = "tukey")

# Two-way interactions
pairs(emmeans(anova_pay, ~ Noise_Type * Control_Condition), adjust = "tukey")
pairs(emmeans(anova_pay, ~ Noise_Type * Noise_Level), adjust = "tukey")
pairs(emmeans(anova_pay, ~ Control_Condition * Noise_Level), adjust = "tukey")

# Three-way interaction
three_way_results <- as.data.frame(
  pairs(
    emmeans(anova_pay, ~ Control_Condition * Noise_Type * Noise_Level),
    adjust = "tukey"
  )
) %>%
  mutate(abs_diff_mean = abs(estimate))

print(three_way_results)

