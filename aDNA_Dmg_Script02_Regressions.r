#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#            Patterns of aDNA Damage Through Time end Environments             #
#                      – lessons from herbarium specimens  -                   #
#                                                                              #
#                                   Script 02                                  #
#                                                                              #
#                        DATA ANALYSIS II - Regressions:                       #
#                                                                              #
#                    - Endogenous Fraction ~ Collection year                   #
#                      - Fragment length ~ Collection year                     #
#                             - 5' C>T  ~ Sample Age                           #
#                  - Damage fraction per site (λ) ~ Sample Age                 #
#------------------------------------------------------------------------------#

# Install required packages, uncomment below if needed
# install.packages(c("dplyr", "ggplot2", "purrr", "gridExtra", "grid", "cowplot"))

# Load libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(gridExtra)
library(grid)
library(cowplot)

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
  mutate(Sample_Age = current_year - Collection_year)

#------------------------------------------------------------------------------#
#          1 - Regression analysis: Endogenous fraction ~ Collection year      #
#------------------------------------------------------------------------------#
# Regression: Endogenous fraction ~ Sample age
# Divided by genera
hordeum_model <- lm(Endogenous_fraction ~ Collection_year,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(Endogenous_fraction ~ Collection_year,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
# Create plot with separate regression lines
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
P1 <- ggplot(d_filtered, aes(x = Collection_year, y = Endogenous_fraction, color = Genus)) +
    geom_point(alpha = 0.5, size = 0.75) +
    geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
    scale_x_continuous(breaks = seq(1800, 2024, by = 20)) +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Collection year",
         y = "% Endogenous DNA",
         title = NULL,
         color = NULL) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom",) +
    # Add R-squared and p-value annotations for both groups separately
    annotate("text", size = 3.5, x = 1800, y = 12,
             label = sprintf("italic(Hordeum):
            italic(R)^2 == %.3f~~italic(p) == %.3g",
                             hordeum_r2, hordeum_p),
             parse = TRUE, hjust = 0) +
    annotate("text", size = 3.5, x = 1800, y = 6,
             label = sprintf("italic(Oryza):
            italic(R)^2 == %.3f~~italic(p) == %.3g",
                             oryza_r2, oryza_p),
             parse = TRUE, hjust = 0)
# Print plot
P1

#------------------------------------------------------------------------------#
#     1.1 - Regression analysis: Endogenous fraction ~ Collection year         #
#                                All Samples                                   #
#------------------------------------------------------------------------------#
# Regression: Endogenous fraction ~ collection year
# For all samples 
model <- lm(Endogenous_fraction ~ Collection_year, data = d_filtered)
# Extract R-squared and p-value
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]
# Define color palette for genera
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
# Plot data
P1.1 <- ggplot(d_filtered, aes(x = Collection_year, y = Endogenous_fraction, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75)+  # Set transparency to 0.5, size to 0.75
  geom_smooth(method = "lm", col = "black") +
  scale_x_continuous(breaks = seq(1800, 2024, by = 20)) +  # Adjust x-axis
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Year",
       y = "% Endogenous DNA",
       title = NULL,
       color = NULL) +  # Remove legend title
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    legend.position = "bottom") +
  # Add R-squared and p-value annotations
  annotate("text", size = 3.5, x = 1970, y = 15,
           label = sprintf("italic(R)^2 == %.3f", r_squared),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 1970, y = 7,
           label = sprintf("italic(p) == %.3g", p_value),
           parse = TRUE, hjust = 0)
# print plot
P1.1

#------------------------------------------------------------------------------#
#         2 - Regression analysis: Fragment length ~ Collection year           #
#                               Divided by Genera                              #
#------------------------------------------------------------------------------#
# Regression: log-mean of fragment length ~ collection year
# Divided by genera
hordeum_model <- lm(Log_Mean ~ Collection_year,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(Log_Mean ~ Collection_year,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
# Create plot with separate regression lines
P2 <- ggplot(d_filtered, aes(x = Collection_year, y = MEDIAN_SIZE, color = Genus)) +
    geom_point(alpha = 0.5, size = 0.75) +
    geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
    scale_y_log10(breaks = c(20, 30, 40, 50, 60, 70),
                  limits = c(19, 70)) +
    scale_x_continuous(breaks = seq(1800, 2024, by = 20)) +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Collection Year",
         y = "Median Fragment Length (log-scale)",
         title = NULL,
         color = NULL) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "bottom") +
    # Add R-squared and p-value annotations for both groups separately
    annotate("text", size = 3.5, x = 1800, y = 21.75,
             label = sprintf("italic(Hordeum):
            italic(R)^2 == %.3f~~italic(p) == %.3g",
                             hordeum_r2, hordeum_p),
             parse = TRUE, hjust = 0) +
    annotate("text", size = 3.5, x = 1800, y = 20,
             label = sprintf("italic(Oryza):
            italic(R)^2 == %.3f~~italic(p) == %.3g",
                             oryza_r2, oryza_p),
             parse = TRUE, hjust = 0)
# print plot
P2

#------------------------------------------------------------------------------#
#         2.1 - Regression analysis: Fragment length ~ Collection year         #
#                                All Samples                                   #
#------------------------------------------------------------------------------#
# Regression: log-mean of fragment length ~ collection year
# For all samples irrespective of Herbarium and genus
model <- lm(Log_Mean ~ Collection_year, data = d_filtered)
# Extract R-squared and p-value
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]
# Define color palette for genera
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
# Plot data
P2.1 <- ggplot(d_filtered, aes(x = Collection_year, y = MEDIAN_SIZE, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75)+  # Set transparency to 0.5, size to 0.75
  geom_smooth(method = "lm", col = "black") +
  scale_y_log10(breaks = c(20, 30, 40, 50, 60, 70),
                limits = c(19, 70)) +  # Set specific y-axis Log-scale
  scale_x_continuous(breaks = seq(1800, 2024, by = 20)) +  # Adjust x-axis
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Year",
       y = "Median Fragment Length (log-scale)",
       title = NULL,
       color = NULL) +  # Remove legend title
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    legend.position = "bottom") +
  # Add R-squared and p-value annotations
  annotate("text", size = 3.5, x = 1970, y = 22.75,
           label = sprintf("italic(R)^2 == %.3f", r_squared),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 1970, y = 20.5,
           label = sprintf("italic(p) == %.3g", p_value),
           parse = TRUE, hjust = 0)
# print plot
P2.1

#------------------------------------------------------------------------------#
#               3 - Regression analysis: 5' C>T  ~ Sample Age                  #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ sample age - by Genera
hordeum_model <- lm(X5P_DMG_POS1 ~ Sample_Age,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_DMG_POS1 ~ Sample_Age,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
P3 <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                      limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "5' C>T Frequencies",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
    legend.position = "bottom") +
  # Add R-squared and p-value annotations for both groups separately
  annotate("text", size = 3.5, x = 130, y = 0.0035,
           label = sprintf("italic(Hordeum): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           hordeum_r2, hordeum_p),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 130, y = 0.000005,
           label = sprintf("italic(Oryza): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           oryza_r2, oryza_p),
           parse = TRUE, hjust = 0)
# print plot
P3

#------------------------------------------------------------------------------#
#               3.1 - Regression analysis: 5' C>T  ~ Sample Age                #
#                                All Samples                                   #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ sample age - all together
model <- lm(X5P_DMG_POS1 ~ Sample_Age, data = d_filtered)
# Extract R-squared and p-value
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]
#Plot 5' damage vs. Sample Age
P3.1 <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_DMG_POS1, color = Genus)) +
    geom_point(alpha = 0.5, size = 0.75) +  # Set transparency to 0.5
    geom_smooth(method = "lm", col = "black") +
    scale_color_manual(values = genus_colors) +  # Set colors for genera
    scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                      limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
    scale_x_continuous(breaks = seq(0, 300, by = 50)) +  # Adjust x-axis
    labs(x = "Sample Age (years)",
         y = "5' C>T frequencies",
         title = NULL,
         color = NULL) +  # Remove legend title
    theme_minimal() +
    theme(panel.grid = element_blank(),  # Remove gridlines
          legend.position = "bottom") +
    # Add R-squared and p-value annotations at the bottom left
    annotate("text", size = 3.5, x = 180, y = 0.0039,
             label = sprintf("italic(R)^2 == %.3f", r_squared),
             parse = TRUE, hjust = 0) +
    annotate("text", size = 3.5, x = 180, y = 0.00014,
             label = sprintf("italic(p) == %.3g", p_value),
             parse = TRUE, hjust = 0)
# print plot
P3.1

#------------------------------------------------------------------------------#
#      4 - Regression analysis: damage fraction per site (λ) and sample age    #
#                           Divided by Genera                                  #
#------------------------------------------------------------------------------#
# Regression: damage fraction per site (λ) ~ sample age
# Divided by Genera
# Format k values in scientific notation in the form 1.2×10^xx
format_scientific <- function(x) {
  # Get the exponent
  exponent <- floor(log10(abs(x)))
  # Get the coefficient (mantissa)
  coefficient <- x / 10^exponent
  # Format as "a.bc×10^xx"
  return(sprintf("%.2f %s 10^%d", coefficient, "%*%", exponent))
}

# Hordeum
hordeum_model <- lm(Lambda ~ Sample_Age,
  data = d_filtered[d_filtered$Genus == "Hordeum", ])
hordeum_r_squared <- summary(hordeum_model)$r.squared
hordeum_p_value <- summary(hordeum_model)$coefficients[2, 4]
hordeum_k <- format_scientific(hordeum_model$coefficients[2])

# Oryza
oryza_model <- lm(Lambda ~ Sample_Age,
  data = d_filtered[d_filtered$Genus == "Oryza", ])
oryza_r_squared <- summary(oryza_model)$r.squared
oryza_p_value <- summary(oryza_model)$coefficients[2, 4]
oryza_k <- format_scientific(oryza_model$coefficients[2])

# Generate Plot
P4 <- ggplot(d_filtered, aes(x = Sample_Age, y = Lambda, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", aes(group = Genus, color = Genus), se = TRUE) +
  scale_color_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "Damage Fraction per Site (λ)",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom" ) +
  annotate(
    "text", size = 3.5, x = 100, y = -0.001,
    label = sprintf("italic(Hordeum):
 italic(R)^2 == %.3f~~italic(p) == %.3g~~k == %s",
     hordeum_r_squared, hordeum_p_value, hordeum_k),
    parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 100, y = -0.012,
    label = sprintf("italic(Oryza):
 italic(R)^2 == %.3f~~italic(p) == %.3g~~k == %s",
     oryza_r_squared, oryza_p_value, oryza_k),
    parse = TRUE, hjust = 0)
# print plot
P4

#------------------------------------------------------------------------------#
#      4.1 - Regression analysis: damage fraction per site (λ) and sample age  #
#                                 All Samples                                  #
#------------------------------------------------------------------------------#
# Regression: damage fraction per site (λ) ~ sample age
# All samples
model <- lm(Lambda ~ Sample_Age, data = d_filtered)
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]
k <- format_scientific(model$coefficients[2])

# Visualization of damage fraction per site (λ) vs. Sample Age
P4.1 <- ggplot(d_filtered, aes(x = Sample_Age, y = Lambda, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", col = "black") +
  scale_color_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "Damage Fraction per Site (λ)",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom") +
  annotate("text", size = 3.5, x = 180, y = 0.018,
           label = sprintf("italic(R)^2 == %.3f", r_squared),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 180, y = 0.008,
           label = sprintf("italic(p) == %.3g", p_value),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 180, y = 0.001,
           label = sprintf("k == %s", k),
           parse = TRUE, hjust = 0)
# print plot
P4.1

#------------------------------------------------------------------------------#
# Print main plots together with labels                                       #
#------------------------------------------------------------------------------#
# 1. Create a separate plot just for the legend
legend_plot <- P1 + theme(legend.position = "bottom")
# 2. Extract just the legend
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
# 3. Remove legends from plots and add labels
P1 <- P1 + theme(legend.position = "none") + labs(tag = "a")
P2 <- P2 + theme(legend.position = "none") + labs(tag = "b")
P4 <- P4 + theme(legend.position = "none") + labs(tag = "c")
P3 <- P3 + theme(legend.position = "none") + labs(tag = "d")

# 4. Create the 2x2 grid with legend at bottom
grid.arrange(
  arrangeGrob(P1, P2, P4, P3, ncol = 2),
  legend,
  heights = c(10, 1),
  ncol = 1
)

#------------------------------------------------------------------------------#
# Print supplementary plots together                                           #
#------------------------------------------------------------------------------#
# 1. Create a separate plot just for the legend
legend_plot <- P1.1 + theme(legend.position = "bottom")
# 2. Extract just the legend
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
# 3. Arrange plots without legends
P1.1 <- P1.1 + theme(legend.position = "none") + labs(tag = "a")
P2.1 <- P2.1 + theme(legend.position = "none") + labs(tag = "b")
P4.1 <- P4.1 + theme(legend.position = "none") + labs(tag = "c")
P3.1 <- P3.1 + theme(legend.position = "none") + labs(tag = "d")

# 4. Create the 2x2 grid with legend at bottom
grid.arrange(
  arrangeGrob(P1.1, P2.1, P4.1, P3.1, ncol = 2),
  legend,
  heights = c(10, 1),
  ncol = 1
)


#------------------------------------------------------------------------------#
# Print detailed linear model results for each damage metric                   #
#------------------------------------------------------------------------------#
# 1. Endogenous fraction ~ Collection Year by Genus
# Hordeum model
hordeum_endo_model <- lm(Endogenous_fraction ~ Collection_year,
                         data = filter(d_filtered, Genus == "Hordeum"))
# Oryza model
oryza_endo_model <- lm(Endogenous_fraction ~ Collection_year,
                       data = filter(d_filtered, Genus == "Oryza"))
# Print detailed results
cat("\n===== ENDOGENOUS FRACTION REGRESSION RESULTS =====\n")
cat("\n----- Hordeum: Endogenous fraction ~ Collection Year -----\n")
print(summary(hordeum_endo_model))
cat("\n----- Oryza: Endogenous fraction ~ Collection Year -----\n")
print(summary(oryza_endo_model))

# 2. Fragment Length ~ Collection Year by Genus
# Hordeum model
hordeum_frag_model <- lm(Log_Mean ~ Collection_year,
                         data = filter(d_filtered, Genus == "Hordeum"))
# Oryza model
oryza_frag_model <- lm(Log_Mean ~ Collection_year,
                       data = filter(d_filtered, Genus == "Oryza"))
# Print detailed results
cat("\n===== FRAGMENT LENGTH REGRESSION RESULTS =====\n")
cat("\n----- Hordeum: Log Mean Fragment Length ~ Collection Year -----\n")
print(summary(hordeum_frag_model))
cat("\n----- Oryza: Log Mean Fragment Length ~ Collection Year -----\n")
print(summary(oryza_frag_model))

# 3. 5' C>T ~ Sample Age by Genus
# Hordeum model
hordeum_ct_model <- lm(X5P_DMG_POS1 ~ Sample_Age,
                       data = filter(d_filtered, Genus == "Hordeum"))
# Oryza model
oryza_ct_model <- lm(X5P_DMG_POS1 ~ Sample_Age,
                     data = filter(d_filtered, Genus == "Oryza"))
# Print detailed results
cat("\n===== 5' C>T DAMAGE REGRESSION RESULTS =====\n")
cat("\n----- Hordeum: 5' C>T Frequencies ~ Sample Age -----\n")
print(summary(hordeum_ct_model))
cat("\n----- Oryza: 5' C>T Frequencies ~ Sample Age -----\n")
print(summary(oryza_ct_model))

# 4. Lambda (Damage fraction per site) ~ Sample Age by Genus
# Hordeum model
hordeum_lambda_model <- lm(Lambda ~ Sample_Age,
                           data = d_filtered[d_filtered$Genus == "Hordeum", ])
# Calculate k value
hordeum_k_value <- hordeum_lambda_model$coefficients[2] * 10^4
# Oryza model
oryza_lambda_model <- lm(Lambda ~ Sample_Age,
                         data = d_filtered[d_filtered$Genus == "Oryza", ])
# Calculate k value
oryza_k_value <- oryza_lambda_model$coefficients[2] * 10^5
# Print detailed results
cat("\n===== DAMAGE FRACTION PER SITE (λ) REGRESSION RESULTS =====\n")
cat("\n----- Hordeum: Damage Fraction per Site (λ) ~ Sample Age -----\n")
print(summary(hordeum_lambda_model))
cat("\nHordeum k value:", sprintf("%.3f × 10^-4", hordeum_k_value), "\n")
cat("\n----- Oryza: Damage Fraction per Site (λ) ~ Sample Age -----\n")
print(summary(oryza_lambda_model))
cat("\nOryza k value:", sprintf("%.3f × 10^-5", oryza_k_value), "\n")

