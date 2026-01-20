#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#            Patterns of aDNA Damage Through Time end Environments             #
#                      – lessons from herbarium specimens  -                   #
#                                                                              #
#                                   Script 04                                  #
#                                                                              #
#                      DATA ANALYSIS IV - COVARIANCE (ANCOVA):                 #
#                                                                              #
#                      - Endogenous fraction ~ Sample Age                      #
#                          - Fragment size ~ Sample Age                        #
#                       - 5' C>T frequencies ~ Sample Age                      #
#                              - Lambda ~ Sample Age                           #
#                                                                              #
#           - Model comparison (interaction vs. no interaction) - ANOVA        #
#                              - Generate Boxplots                             #
#                                                                              #
#------------------------------------------------------------------------------#

# Install required packages, uncomment below if needed
# install.packages(c("dplyr", "ggplot2", "ggpubr", "car"))

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(car)

# Read data matrix in R
d <- read.delim("aDNA_damage_screening_MAIN.txt",
                row.names = 1, stringsAsFactors = TRUE)
# Filter ( minimum 5000 merged reads, remove post-2010,
# remove samples deviating from exponential decline of fragment length and 5'Damage)
# Define the filtering criteria
filter_condition <- d$NO_MERGED_READS >= 5000 &
  d$Collection_year < 2010 &
  d$Lambda_R_squared > 0.95 &
  d$Lambda_p_value < 0.05 &
  d$X5P_DMG_R_squared > 0.5 &
  d$X5P_DMG_p_value < 0.05

# Filter the data to keep samples meeting all criteria
d_filtered <- d[filter_condition, ]
# Store samples that did not meet the criteria
d_removed <- d[!filter_condition, ]

# Filter out rows with non-finite or zero values in
# MEDIAN_SIZE, Log_Mean or Collection_year
# Calculate age of the sample
current_year <- 2025
d_filtered <- d_filtered %>%
  filter(!is.na(MEDIAN_SIZE) & MEDIAN_SIZE > 0,
         !is.na(Log_Mean),
         !is.na(Collection_year)) %>%
  mutate(Sample_Age = current_year - Collection_year,
         # Create corrected 5' C>T value by subtracting non-deamination background
         X5P_DMG_POS1_Corrected = X5P_DMG_POS1 - X5P_other_freq,
         X3P_DMG_POS1_Corrected = X3P_DMG_POS1 - X3P_other_freq)

# Set color palette for genera
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")

#------------------------------------------------------------------------------#
# 1. Endogenous Fraction - Model Comparison                                    #
#------------------------------------------------------------------------------#
# Full model with interaction
endo_full <- lm(Endogenous_fraction ~ Sample_Age * Genus, data = d_filtered)
summary_endo_full <- summary(endo_full)
anova_endo_full <- anova(endo_full)

# Reduced model without interaction
endo_reduced <- lm(Endogenous_fraction ~ Sample_Age + Genus, data = d_filtered)
summary_endo_reduced <- summary(endo_reduced)

# Compare models using ANOVA
endo_comparison <- anova(endo_reduced, endo_full)
print("Model comparison for Endogenous Fraction:")
print(endo_comparison)

# Extract interaction p-value for model selection
endo_interaction_p <- anova_endo_full["Sample_Age:Genus", "Pr(>F)"]

# Select final model based on interaction significance
endo_final <- if(endo_interaction_p < 0.05) endo_full else endo_reduced
endo_final_summary <- summary(endo_final)
endo_final_type <- if(endo_interaction_p < 0.05) "with interaction" else "without interaction"

# Extract p-values from the final model for annotation
if(endo_interaction_p < 0.05) {
  # Use p-values from full model if interaction is significant
  endo_age_p <- summary_endo_full$coefficients["Sample_Age", "Pr(>|t|)"]
  endo_genus_p <- summary_endo_full$coefficients["GenusOryza", "Pr(>|t|)"]
  endo_int_p <- endo_interaction_p
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          endo_age_p, endo_genus_p, endo_int_p)
} else {
  # Use p-values from reduced model if interaction is not significant
  endo_age_p <- summary_endo_reduced$coefficients["Sample_Age", "Pr(>|t|)"]
  endo_genus_p <- summary_endo_reduced$coefficients["GenusOryza", "Pr(>|t|)"]
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          endo_age_p, endo_genus_p, endo_interaction_p)
}

# Visualization - Boxplot
endo_box <- ggplot(d_filtered, aes(x = Genus, y = Endogenous_fraction, fill = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_x_discrete(labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "% Endogenous DNA",
       title = NULL,
       caption = caption_text) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#------------------------------------------------------------------------------#
# 2. Fragment Length (Log_Mean) - Model Comparison                             #
#------------------------------------------------------------------------------#
# Full model with interaction
frag_full <- lm(Log_Mean ~ Sample_Age * Genus, data = d_filtered)
summary_frag_full <- summary(frag_full)
anova_frag_full <- anova(frag_full)

# Reduced model without interaction
frag_reduced <- lm(Log_Mean ~ Sample_Age + Genus, data = d_filtered)
summary_frag_reduced <- summary(frag_reduced)

# Compare models using ANOVA
frag_comparison <- anova(frag_reduced, frag_full)
print("Model comparison for Fragment Length:")
print(frag_comparison)

# Extract interaction p-value for model selection
frag_interaction_p <- anova_frag_full["Sample_Age:Genus", "Pr(>F)"]

# Select final model based on interaction significance
frag_final <- if(frag_interaction_p < 0.05) frag_full else frag_reduced
frag_final_summary <- summary(frag_final)
frag_final_type <- if(frag_interaction_p < 0.05) "with interaction" else "without interaction"

# Extract p-values from the final model for annotation
if(frag_interaction_p < 0.05) {
  # Use p-values from full model if interaction is significant
  frag_age_p <- summary_frag_full$coefficients["Sample_Age", "Pr(>|t|)"]
  frag_genus_p <- summary_frag_full$coefficients["GenusOryza", "Pr(>|t|)"]
  frag_int_p <- frag_interaction_p
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          frag_age_p, frag_genus_p, frag_int_p)
} else {
  # Use p-values from reduced model if interaction is not significant
  frag_age_p <- summary_frag_reduced$coefficients["Sample_Age", "Pr(>|t|)"]
  frag_genus_p <- summary_frag_reduced$coefficients["GenusOryza", "Pr(>|t|)"]
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          frag_age_p, frag_genus_p, frag_interaction_p)
}

# Visualization - Boxplot 
frag_box <- ggplot(d_filtered, aes(x = Genus, y = MEDIAN_SIZE, fill = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 21, outlier.size = 1.5) +
  scale_y_log10(breaks = c(30, 40, 50, 60, 70),
                limits = c(30, 70)) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_x_discrete(labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "Median Fragment Length (log-scaled)",
       title = NULL,
       caption = caption_text) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#------------------------------------------------------------------------------#
# 3. 5' C>T Damage (Corrected) - Model Comparison                              #
#------------------------------------------------------------------------------#
# Full model with interaction
damage_full <- lm(X5P_DMG_POS1 ~ Sample_Age * Genus, data = d_filtered)
summary_damage_full <- summary(damage_full)
anova_damage_full <- anova(damage_full)

# Reduced model without interaction
damage_reduced <- lm(X5P_DMG_POS1 ~ Sample_Age + Genus, data = d_filtered)
summary_damage_reduced <- summary(damage_reduced)

# Compare models using ANOVA
damage_comparison <- anova(damage_reduced, damage_full)
print("Model comparison for 5' C>T Damage:")
print(damage_comparison)

# Extract interaction p-value for model selection
damage_interaction_p <- anova_damage_full["Sample_Age:Genus", "Pr(>F)"]

# Select final model based on interaction significance
damage_final <- if(damage_interaction_p < 0.05) damage_full else damage_reduced
damage_final_summary <- summary(damage_final)
damage_final_type <- if(damage_interaction_p < 0.05) "with interaction" else "without interaction"

# Extract p-values from the final model for annotation
if(damage_interaction_p < 0.05) {
  # Use p-values from full model if interaction is significant
  damage_age_p <- summary_damage_full$coefficients["Sample_Age", "Pr(>|t|)"]
  damage_genus_p <- summary_damage_full$coefficients["GenusOryza", "Pr(>|t|)"]
  damage_int_p <- damage_interaction_p
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          damage_age_p, damage_genus_p, damage_int_p)
} else {
  # Use p-values from reduced model if interaction is not significant
  damage_age_p <- summary_damage_reduced$coefficients["Sample_Age", "Pr(>|t|)"]
  damage_genus_p <- summary_damage_reduced$coefficients["GenusOryza", "Pr(>|t|)"]
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          damage_age_p, damage_genus_p, damage_interaction_p)
}

# Visualization - Boxplot
damage_box <- ggplot(d_filtered, aes(x = Genus, y = X5P_DMG_POS1, fill = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_x_discrete(labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "5' C>T Frequencies",
       title = NULL,
       caption = caption_text) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#------------------------------------------------------------------------------#
# 4. Lambda - Model Comparison                                                 #
#------------------------------------------------------------------------------#
# Full model with interaction
lambda_full <- lm(Lambda ~ Sample_Age * Genus, data = d_filtered)
summary_lambda_full <- summary(lambda_full)
anova_lambda_full <- anova(lambda_full)

# Reduced model without interaction
lambda_reduced <- lm(Lambda ~ Sample_Age + Genus, data = d_filtered)
summary_lambda_reduced <- summary(lambda_reduced)

# Compare models using ANOVA
lambda_comparison <- anova(lambda_reduced, lambda_full)
print("Model comparison for Lambda:")
print(lambda_comparison)

# Extract interaction p-value for model selection
lambda_interaction_p <- anova_lambda_full["Sample_Age:Genus", "Pr(>F)"]

# Select final model based on interaction significance
lambda_final <- if(lambda_interaction_p < 0.05) lambda_full else lambda_reduced
lambda_final_summary <- summary(lambda_final)
lambda_final_type <- if(lambda_interaction_p < 0.05) "with interaction" else "without interaction"

# Extract p-values from the final model for annotation
if(lambda_interaction_p < 0.05) {
  # Use p-values from full model if interaction is significant
  lambda_age_p <- summary_lambda_full$coefficients["Sample_Age", "Pr(>|t|)"]
  lambda_genus_p <- summary_lambda_full$coefficients["GenusOryza", "Pr(>|t|)"]
  lambda_int_p <- lambda_interaction_p
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          lambda_age_p, lambda_genus_p, lambda_int_p)
} else {
  # Use p-values from reduced model if interaction is not significant
  lambda_age_p <- summary_lambda_reduced$coefficients["Sample_Age", "Pr(>|t|)"]
  lambda_genus_p <- summary_lambda_reduced$coefficients["GenusOryza", "Pr(>|t|)"]
  caption_text <- sprintf("Age effect: p = %.3g, Genus effect: p = %.3g, Interaction: p = %.3g",
                          lambda_age_p, lambda_genus_p, lambda_interaction_p)
}

# Visualization - Boxplot
lambda_box <- ggplot(d_filtered, aes(x = Genus, y = Lambda, fill = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_x_discrete(labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "Damage Fraction per Site (λ)",
       title = NULL,
       caption = caption_text) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#------------------------------------------------------------------------------#
# Combine boxplots in a single figure                                          #
#------------------------------------------------------------------------------#
combined_box <- ggpubr::ggarrange(
  endo_box, frag_box, lambda_box, damage_box,
  ncol = 2, nrow = 2, 
  common.legend = TRUE,
  legend = "bottom"
)

# Print the combined figure to screen
print(combined_box)

#------------------------------------------------------------------------------#
# Print full summary table of all models & model comparisons                   #
#------------------------------------------------------------------------------#

# Endogenoun fraction
summary(endo_full)
summary(endo_reduced)
print(endo_comparison)

# Fragment size
summary(frag_full)
summary(frag_reduced)
print(frag_comparison)

# 5' C>T
summary(damage_full)
summary(damage_reduced)
print(damage_comparison)

# lambda
summary(lambda_full)
summary(lambda_reduced)
print(lambda_comparison)

#------------------------------------------------------------------------------#
# Print a summary table of all model comparisons                              #
#------------------------------------------------------------------------------#
# Create a data frame with model comparison results
model_summary <- data.frame(
  Parameter = c("Endogenous Fraction", "Fragment Length", "5' C>T Damage", "Lambda"),
  Interaction_P = c(endo_interaction_p, frag_interaction_p, damage_interaction_p, lambda_interaction_p),
  Model_Comparison_P = c(endo_comparison$`Pr(>F)`[2], 
                         frag_comparison$`Pr(>F)`[2],
                         damage_comparison$`Pr(>F)`[2],
                         lambda_comparison$`Pr(>F)`[2]),
  Final_Model = c(endo_final_type, frag_final_type, damage_final_type, lambda_final_type),
  R_squared_full = c(summary(endo_full)$r.squared,
                     summary(frag_full)$r.squared,
                     summary(damage_full)$r.squared,
                     summary(lambda_full)$r.squared),
  R_squared_reduced = c(summary(endo_reduced)$r.squared,
                       summary(frag_reduced)$r.squared,
                       summary(damage_reduced)$r.squared,
                       summary(lambda_reduced)$r.squared),
  R_squared_final = c(summary(endo_final)$r.squared,
                     summary(frag_final)$r.squared,
                     summary(damage_final)$r.squared,
                     summary(lambda_final)$r.squared)
)

# Format p-values and R-squared values
model_summary$Interaction_P <- formatC(model_summary$Interaction_P, format = "e", digits = 3)
model_summary$Model_Comparison_P <- formatC(model_summary$Model_Comparison_P, format = "e", digits = 3)
model_summary$R_squared_full <- round(model_summary$R_squared_full, 4)
model_summary$R_squared_reduced <- round(model_summary$R_squared_reduced, 4)
model_summary$R_squared_final <- round(model_summary$R_squared_final, 4)

# Print summary table
print("Summary of Model Comparisons:")
print(model_summary)

# Write summary table to CSV
write.csv(model_summary, "model_comparison_summary.csv", row.names = FALSE)
