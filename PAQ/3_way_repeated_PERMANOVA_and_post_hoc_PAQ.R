# Load required libraries
library(lme4)
library(carData)
library(car)
library(emmeans)
library(dplyr)
library(multcomp)
library(tidyr)
library(vegan)
library(pairwiseAdonis)

# --------------------------
# Section 2: Data preparation
# --------------------------

# Read data, replacing empty strings with NA
paq_data <- read.csv("PAQ_results.csv", stringsAsFactors = TRUE, na.strings = c("", "NA"))
stimuli_data <- read.csv("Stimuli_information.csv", stringsAsFactors = TRUE)

# Rename columns for clarity
paq_data <- paq_data %>% 
  rename_with(~ gsub("Stimulus\\.(\\d+)\\.\\.(X|Y)\\.$", "Stimuli_\\1_\\2", .))

# Add a column for participant numbers (if not already included)
if (!"Participants" %in% colnames(paq_data)) {
  paq_data <- paq_data %>% mutate(Participants = row_number())
}

# Convert data to long format
paq_scores <- paq_data %>%
  pivot_longer(
    cols = starts_with("Stimuli_"),   # Focus on stimulus columns
    names_to = "Stimulus",            # Name of transformed columns
    values_to = "Score"               # Name of associated values
  ) %>%
  mutate(
    Stimulus_Num = gsub("Stimuli_(\\d+)_.*", "\\1", Stimulus),  # Extract stimulus number
    Dimension = gsub(".*_(X|Y)", "\\1", Stimulus)              # Identify X or Y
  ) %>%
  pivot_wider(
    names_from = Dimension, 
    values_from = Score      # Spread scores into X and Y columns
  ) %>%
  dplyr::select(Participants, Stimulus_Num, X, Y) # Keep only relevant columns

# Group by participant and stimulus, then take the first non-NA value for each dimension
paq_scores <- paq_scores %>%
  group_by(Participants, Stimulus_Num) %>%
  summarise(
    X = first(na.omit(X)),  # Use the first non-NA value for X
    Y = first(na.omit(Y))   # Use the first non-NA value for Y
  ) %>%
  ungroup()

# Add the prefix "Stimuli_" to Stimulus_Num to match the "Stimuli" column in stimuli_data
paq_scores$Stimulus_Num <- paste0("Stimuli_", paq_scores$Stimulus_Num)

# Merge data without creating additional columns for Noise_Type, Control_Condition, and Noise_Level
paq_scores <- paq_scores %>%
  left_join(stimuli_data, by = c("Stimulus_Num" = "Stimuli"))

# Check the result
head(paq_scores)

# --------------------------
# Section 3: Non-parametric repeated-measures MANOVA (PERMANOVA)
# --------------------------

# Ensure factors are correctly specified as factors
paq_scores$Noise_Type <- as.factor(paq_scores$Noise_Type)
paq_scores$Control_Condition <- as.factor(paq_scores$Control_Condition)
paq_scores$Noise_Level <- as.factor(paq_scores$Noise_Level)
paq_scores$Participants <- as.factor(paq_scores$Participants)

# Remove rows with missing values in X or Y
paq_scores_clean <- paq_scores %>%
  dplyr::filter(!is.na(X) & !is.na(Y))

# Create a data matrix with X and Y as dependent variables
score_matrix <- cbind(paq_scores_clean$X, paq_scores_clean$Y)

# Create the distance matrix (Euclidean distance)
distance_matrix <- dist(score_matrix)

# Run PERMANOVA tests for different configurations
permanova_results <- list()

# Create a list of formulas for PERMANOVA models
models <- list(
  Noise_Type = "Noise_Type + Participants",
  Control_Condition = "Control_Condition + Participants",
  Noise_Level = "Noise_Level + Participants",
  Noise_Type_Control_Condition = "Noise_Type * Control_Condition + Participants",
  Control_Condition_Noise_Level = "Control_Condition * Noise_Level + Participants",
  Noise_Type_Noise_Level = "Noise_Type * Noise_Level + Participants",
  Noise_Type_Control_Condition_Noise_Level = "Noise_Type * Control_Condition * Noise_Level + Participants"
)

# Initialize an empty list to store results
permanova_results <- list()

# Run PERMANOVA for each formula in the list
for (model_name in names(models)) {
  formula <- as.formula(paste("distance_matrix ~", models[[model_name]]))
  permanova_results[[model_name]] <- adonis2(formula, data = paq_scores_clean, permutations = 999)
}

# Extract results into a table
results_table <- lapply(permanova_results, function(result) {
  data.frame(
    Df = result$aov.tab$Df,              # Degrees of freedom
    SumOfSqs = result$aov.tab$`Sum Sq`,  # Sum of squares
    R2 = result$aov.tab$R2,              # R-squared
    F = result$aov.tab$`F value`,        # F-statistic
    Pr_F = result$aov.tab$`Pr(>F)`       # p-value
  )
})

# Combine results into a single table
final_results <- bind_rows(results_table, .id = "Factor")

# Display the final table
print(final_results)

# Create an empty table to store results
results_table <- data.frame(Configuration = character(),
                            F_value = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Update interaction names
two_way_interactions <- c("Noise_Type_Control_Condition", "Noise_Type_Noise_Level", "Control_Condition_Noise_Level")
three_way_interaction <- "Noise_Type_Control_Condition_Noise_Level"

# Simple configurations
simple_configurations <- c("Noise_Type", "Control_Condition", "Noise_Level")

# Initialize an empty table
results_table <- data.frame(Configuration = character(),
                            F_value = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Add simple configurations
for (config in simple_configurations) {
  F_value <- permanova_results[[config]]$F[1]
  p_value <- permanova_results[[config]]$`Pr(>F)`[1]
  if (!config %in% results_table$Configuration) {
    results_table <- rbind(results_table, data.frame(Configuration = config,
                                                     F_value = F_value,
                                                     p_value = p_value))
  }
}

# Add two-way interactions
for (interaction in two_way_interactions) {
  F_value <- permanova_results[[interaction]]$F[1]
  p_value <- permanova_results[[interaction]]$`Pr(>F)`[1]
  if (!interaction %in% results_table$Configuration) {
    results_table <- rbind(results_table, data.frame(Configuration = interaction,
                                                     F_value = F_value,
                                                     p_value = p_value))
  }
}

# Add the three-way interaction
F_value <- permanova_results[[three_way_interaction]]$F[1]
p_value <- permanova_results[[three_way_interaction]]$`Pr(>F)`[1]
if (!three_way_interaction %in% results_table$Configuration) {
  results_table <- rbind(results_table, data.frame(Configuration = three_way_interaction,
                                                   F_value = F_value,
                                                   p_value = p_value))
}

# Filter rows with NA F or p values
results_table_cleaned <- results_table[!is.na(results_table$F_value) & !is.na(results_table$p_value), ]

# Apply FDR correction to p-values
results_table_cleaned$p_value_fdr <- p.adjust(results_table_cleaned$p_value, method = "fdr")

# Display results with corrected p-values
print(results_table_cleaned)

# Display the cleaned results table
print(results_table_cleaned)

# Post-hoc analysis

# Function to calculate means and pairwise differences
calculate_pairwise_differences <- function(data, factor_column, pairwise_results) {
  # Compute means per factor level
  means <- aggregate(data ~ factor_column, data = data.frame(data, factor_column), FUN = mean)
  colnames(means)[1] <- "Factor_Level"
  
  # Create vectors of means and factor levels
  means_vector <- means$data
  factor_levels <- means$Factor_Level
  
  # Add column for mean differences
  pairwise_results$mean_differences <- NA
  
  # Compute mean differences for each pair
  for (i in 1:nrow(pairwise_results)) {
    group_1 <- strsplit(pairwise_results$pairs[i], " vs ")[[1]][1]
    group_2 <- strsplit(pairwise_results$pairs[i], " vs ")[[1]][2]
    
    mean_1 <- means_vector[which(factor_levels == group_1)]
    mean_2 <- means_vector[which(factor_levels == group_2)]
    
    pairwise_results$mean_differences[i] <- (mean_1 - mean_2)
  }
  
  return(pairwise_results)
}

# Function to filter results by adjusted p-value and mean difference
filter_results <- function(pairwise_results) {
  filtered_results <- pairwise_results[
    pairwise_results$p.adjusted < 0.05 & abs(pairwise_results$mean_differences) > 0.25, 
  ]
  return(filtered_results)
}

# 1. Pairwise comparison for Noise_Type
factors_noise <- paq_scores_clean$Noise_Type
pairwise_results_noise <- pairwise.adonis(distance_matrix, factors = factors_noise, perm = 999)
pairwise_results_noise$p.values$`p.value` <- p.adjust(pairwise_results_noise$p.values$`p.value`, method = "BH")
pairwise_results_noise <- calculate_pairwise_differences(paq_scores_clean$X, factors_noise, pairwise_results_noise)
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_noise <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_noise, pairwise_results_noise)
pairwise_results_noise <- filter_results(pairwise_results_noise)
print(pairwise_results_noise)

# 2. Pairwise comparison for Control_Condition
factors_control <- paq_scores_clean$Control_Condition
pairwise_results_control <- pairwise.adonis(distance_matrix, factors = factors_control, perm = 999)
pairwise_results_control$p.values$`p.value` <- p.adjust(pairwise_results_control$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_control <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_control, pairwise_results_control)
pairwise_results_control <- filter_results(pairwise_results_control)
print(pairwise_results_control)

# 3. Pairwise comparison for Noise_Level
factors_level <- paq_scores_clean$Noise_Level
pairwise_results_level <- pairwise.adonis(distance_matrix, factors = factors_level, perm = 999)
pairwise_results_level$p.values$`p.value` <- p.adjust(pairwise_results_level$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_level <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_level, pairwise_results_level)
pairwise_results_level <- filter_results(pairwise_results_level)
print(pairwise_results_level)

# 4. Pairwise comparison for Noise_Type * Control_Condition
factors_2 <- with(paq_scores_clean, interaction(Noise_Type, Control_Condition))
pairwise_results_2 <- pairwise.adonis(distance_matrix, factors = factors_2, perm = 999)
pairwise_results_2$p.values$`p.value` <- p.adjust(pairwise_results_2$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_2 <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_2, pairwise_results_2)
pairwise_results_2 <- filter_results(pairwise_results_2)
print(pairwise_results_2)

# 5. Pairwise comparison for Noise_Type * Noise_Level
factors_2 <- with(paq_scores_clean, interaction(Noise_Type, Noise_Level))
pairwise_results_2 <- pairwise.adonis(distance_matrix, factors = factors_2, perm = 999)
pairwise_results_2$p.values$`p.value` <- p.adjust(pairwise_results_2$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_2 <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_2, pairwise_results_2)
pairwise_results_2 <- filter_results(pairwise_results_2)
print(pairwise_results_2)

# 6. Pairwise comparison for Control_Condition * Noise_Level
factors_2 <- with(paq_scores_clean, interaction(Control_Condition, Noise_Level))
pairwise_results_2 <- pairwise.adonis(distance_matrix, factors = factors_2, perm = 999)
pairwise_results_2$p.values$`p.value` <- p.adjust(pairwise_results_2$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_2 <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_2, pairwise_results_2)
pairwise_results_2 <- filter_results(pairwise_results_2)
print(pairwise_results_2)

# 7. Pairwise comparison for Noise_Type * Control_Condition * Noise_Level
factors_3 <- with(paq_scores_clean, interaction(Noise_Type, Control_Condition, Noise_Level))
pairwise_results_3 <- pairwise.adonis(distance_matrix, factors = factors_3, perm = 999)
pairwise_results_3$p.values$`p.value` <- p.adjust(pairwise_results_3$p.values$`p.value`, method = "BH")
paq_scores_clean$combined_sd <- apply(paq_scores_clean[, c("X", "Y")], 1, sd)
pairwise_results_3 <- calculate_pairwise_differences(paq_scores_clean$combined_sd, factors_3, pairwise_results_3)
pairwise_results_3 <- filter_results(pairwise_results_3)
print(pairwise_results_3)
