#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#                 aDNA Damage through time -  Script 07                        #
#                     SUPPLEMENTARY DEAMINATION ANALYSES                       #
#------------------------------------------------------------------------------#
# Load libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(terra)
library(viridis)
library(gridExtra)
library(vegan)
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
  mutate(Sample_Age = current_year - Collection_year,
         # Create corrected 5' C>T value by subtracting non-deamination background
         X5P_DMG_POS1_Corrected = X5P_DMG_POS1 - X5P_other_freq,
         X3P_DMG_POS1_Corrected = X3P_DMG_POS1 - X3P_other_freq)

# Define genus colors
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")


#------------------------------------------------------------------------------#
#             1. Calculate Divergence from Reference Genome                    #
#   Divergence = Sum of 5' and 3' "other" (non-deamination) frequencies        #
#                                                                              #
#------------------------------------------------------------------------------#

# Calculate baseline substitution rates = evolutionary distance from reference
d_filtered$Divergence_from_Reference <- as.numeric(d_filtered$X5P_other_freq) + 
  as.numeric(d_filtered$X3P_other_freq)

# Summary Statistics by Species
divergence_summary <- d_filtered %>%
  group_by(Species, Genus) %>%
  summarise(
    n_samples = n(),
    Mean_Divergence_pct = mean(Divergence_from_Reference, na.rm = TRUE) * 100,
    SD_Divergence_pct = sd(Divergence_from_Reference, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(Genus, Species)

cat("\n=== Divergence from Reference Genome by Species ===\n")
print(divergence_summary, n = Inf)

#------------------------------------------------------------------------------#
#                              Boxplot                                         #
#------------------------------------------------------------------------------#

Suppl_P1_Divergence_from_Ref <- ggplot(d_filtered, aes(x = Species, y = Divergence_from_Reference, fill = Genus)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Remove default outliers
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +  # Add all data points
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +  # Show as percentages
  labs(x = NULL,
       y = "Divergence from Reference Genome (%)",
       title = NULL,
       fill = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(Suppl_P1_Divergence_from_Ref)

#------------------------------------------------------------------------------#
#                       2. Supplementary Regression analysis:                  #
#                                                                              #
#                             - 5' C>T ~ 3' G>A (RAW)                          #
#             - 5' C>T ~ 3' G>A (corrected by baseline substitution)           #
#                                                                              #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ 3pGA - RAW
hordeum_model <- lm(X5P_DMG_POS1 ~ X3P_DMG_POS1,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_DMG_POS1 ~ X3P_DMG_POS1,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
Suppl_P2_CT_GA_Raw <- ggplot(d_filtered, aes(x = X3P_DMG_POS1, y = X5P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
  scale_x_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     limits = c(0, 0.055)) +  # Adjust x-axis with proper increments
  
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "3' G>A Frequencies",
       y = "5' C>T Frequencies",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") +
  # Add R-squared and p-value annotations for both groups separately
  annotate("text", size = 3.5, x = 0.030, y = 0.0035,
           label = sprintf("italic(Hordeum): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           hordeum_r2, hordeum_p),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 0.030, y = 0.000005,
           label = sprintf("italic(Oryza): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           oryza_r2, oryza_p),
           parse = TRUE, hjust = 0)
# print plot
Suppl_P2_CT_GA_Raw

# Regression: 5pCT ~ 3pGA - Corrected by baseline substitutions
hordeum_model <- lm(X5P_DMG_POS1_Corrected ~ X3P_DMG_POS1_Corrected,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_DMG_POS1_Corrected ~ X3P_DMG_POS1_Corrected,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
Suppl_P2_CT_GA_Corrected <- ggplot(d_filtered, aes(x = X3P_DMG_POS1_Corrected, y = X5P_DMG_POS1_Corrected, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
  scale_x_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     limits = c(0, 0.055)) +  # Adjust x-axis with proper increments
  
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "3' G>A Frequencies (Corrected)",
       y = "5' C>T Frequencies (Corrected)",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") +
  # Add R-squared and p-value annotations for both groups separately
  annotate("text", size = 3.5, x = 0.030, y = 0.0035,
           label = sprintf("italic(Hordeum): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           hordeum_r2, hordeum_p),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 0.030, y = 0.000005,
           label = sprintf("italic(Oryza): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           oryza_r2, oryza_p),
           parse = TRUE, hjust = 0)
# print plot
Suppl_P2_CT_GA_Corrected

#------------------------------------------------------------------------------#
#               3 - Regression analysis: 5' Others  ~ Sample Age               #
#             (to check patterns are specific to deamination rates)            #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ sample age - by Genera
hordeum_model <- lm(X5P_other_freq ~ Sample_Age,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_other_freq ~ Sample_Age,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
Suppl_P3_5p_Others_vs_age <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_other_freq, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "5' Other Frequencies",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") +
  # Add R-squared and p-value annotations for both groups separately
  annotate("text", size = 3.5, x = 130, y = 0.0012,
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
Suppl_P3_5p_Others_vs_age

#------------------------------------------------------------------------------#
#           4 - Regression analysis: 5' C>T (Divergence-Corrected) ~ Age       #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 5' C>T  ~ sample age - by Genera
hordeum_model <- lm(X5P_DMG_POS1_Corrected ~ Sample_Age,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_DMG_POS1_Corrected ~ Sample_Age,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
Suppl_P4_5CT_Corrected_Age <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "5' C>T Frequencies (Corrected)",
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
Suppl_P4_5CT_Corrected_Age

#------------------------------------------------------------------------------#
#           4.1 5' C>T Damage (Corrected) - Model Comparison                   #
#------------------------------------------------------------------------------#
# Full model with interaction
damage_full <- lm(X5P_DMG_POS1_Corrected ~ Sample_Age * Genus, data = d_filtered)
summary_damage_full <- summary(damage_full)
anova_damage_full <- anova(damage_full)

# Reduced model without interaction
damage_reduced <- lm(X5P_DMG_POS1_Corrected ~ Sample_Age + Genus, data = d_filtered)
summary_damage_reduced <- summary(damage_reduced)

# Compare models using ANOVA
damage_comparison <- anova(damage_reduced, damage_full)
print("Model comparison for 5' C>T Damage (Corrected):")
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
Suppl_P5_5CT_Corrected <- ggplot(d_filtered, aes(x = Genus, y = X5P_DMG_POS1_Corrected, fill = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                               expression(italic("Oryza")))) +
  scale_x_discrete(labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "5' C>T Frequencies (Corrected)",
       title = NULL,
       caption = caption_text) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
# print plot
Suppl_P5_5CT_Corrected

#------------------------------------------------------------------------------#
#           4.2 5' C>T Damage (Corrected) - Variance Partitioning              #
#------------------------------------------------------------------------------#

# Clean coordinates
clean_coords <- d_filtered %>%
  filter(!is.na(Latitude_DD), !is.na(Longitude_DD)) %>%
  filter(Latitude_DD >= -90, Latitude_DD <= 90,
         Longitude_DD >= -180, Longitude_DD <= 180)

# Print sample info
print(paste("Number of samples:", nrow(clean_coords)))
print(paste("Samples by genus:", table(clean_coords$Genus)))

# Create points for extraction
sample_points <- vect(clean_coords, geom=c("Longitude_DD", "Latitude_DD"), 
                      crs="epsg:4326")

# Extract bioclimatic values
bio1 <- rast("chelsa_data/CHELSA_bio1_1981-2010_V.2.1.tif")    # Annual Mean Temperature
bio4 <- rast("chelsa_data/CHELSA_bio4_1981-2010_V.2.1.tif")    # Temperature Seasonality
bio12 <- rast("chelsa_data/CHELSA_bio12_1981-2010_V.2.1.tif")  # Annual Precipitation
bio15 <- rast("chelsa_data/CHELSA_bio15_1981-2010_V.2.1.tif")  # Precipitation Seasonality

climate_values <- data.frame(
  temp_annual = extract(bio1, sample_points)[,2],
  temp_seasonality = extract(bio4, sample_points)[,2],
  precip_annual = extract(bio12, sample_points)[,2],
  precip_seasonality = extract(bio15, sample_points)[,2]
)

# Function to load monthly climate data
load_monthly_chelsa <- function(base_path, variable) {
  months <- sprintf("%02d", 1:12)
  paths <- paste0(base_path, "/monthly/", variable, "/CHELSA_", 
                  variable, "_", months, "_1981-2010_V.2.1.tif")
  raster_list <- lapply(paths, rast)
  return(raster_list)
}

# Load monthly temperature (tas) and precipitation (pr) data
monthly_temp <- load_monthly_chelsa("chelsa_data", "tas")
names(monthly_temp) <- paste0("temp_", sprintf("%02d", 1:12))

monthly_precip <- load_monthly_chelsa("chelsa_data", "pr")
names(monthly_precip) <- paste0("precip_", sprintf("%02d", 1:12))

# Extract monthly values
extract_monthly_values <- function(raster_list, points) {
  values_list <- lapply(raster_list, function(x) extract(x, points)[,2])
  values_df <- as.data.frame(do.call(cbind, values_list))
  names(values_df) <- names(raster_list)
  return(values_df)
}

monthly_temp_values <- extract_monthly_values(monthly_temp, sample_points)
monthly_precip_values <- extract_monthly_values(monthly_precip, sample_points)

# Combine all data
data_with_climate <- bind_cols(
  clean_coords,
  climate_values,
  monthly_temp_values,
  monthly_precip_values
)

# Add climate data for collection month
data_analysis <- data_with_climate %>%
  mutate(
    # Get climate data for collection month
    Collection_Temp = case_when(
      Collection_month == 1 ~ temp_01,
      Collection_month == 2 ~ temp_02,
      Collection_month == 3 ~ temp_03,
      Collection_month == 4 ~ temp_04,
      Collection_month == 5 ~ temp_05,
      Collection_month == 6 ~ temp_06,
      Collection_month == 7 ~ temp_07,
      Collection_month == 8 ~ temp_08,
      Collection_month == 9 ~ temp_09,
      Collection_month == 10 ~ temp_10,
      Collection_month == 11 ~ temp_11,
      Collection_month == 12 ~ temp_12
    ),
    Collection_Precip = case_when(
      Collection_month == 1 ~ precip_01,
      Collection_month == 2 ~ precip_02,
      Collection_month == 3 ~ precip_03,
      Collection_month == 4 ~ precip_04,
      Collection_month == 5 ~ precip_05,
      Collection_month == 6 ~ precip_06,
      Collection_month == 7 ~ precip_07,
      Collection_month == 8 ~ precip_08,
      Collection_month == 9 ~ precip_09,
      Collection_month == 10 ~ precip_10,
      Collection_month == 11 ~ precip_11,
      Collection_month == 12 ~ precip_12
    )
  )

# Convert categorical variables to factors
data_analysis$Genus <- as.factor(data_analysis$Genus)
data_analysis$Herbarium <- as.factor(data_analysis$Herbarium)

# Main Variance Partitioning Function
run_main_varpart_analysis <- function(data, response_var, model_type) {
  
  cat("\n", rep("=", 70), "\n")
  cat("ANALYSIS:", response_var, "-", toupper(model_type), "CLIMATE\n")
  cat(rep("=", 70), "\n")
  
  # Filter data and prepare variables based on model type
  if (model_type == "collection") {
    complete_data <- data %>%
      filter(!is.na(!!sym(response_var)),
             !is.na(Sample_Age),
             !is.na(Collection_Temp),
             !is.na(Collection_Precip),
             !is.na(Genus))
    
    X1 <- scale(complete_data$Collection_Temp)
    X2 <- scale(complete_data$Collection_Precip)
    
  } else if (model_type == "annual") {
    complete_data <- data %>%
      filter(!is.na(!!sym(response_var)),
             !is.na(Sample_Age),
             !is.na(temp_annual),
             !is.na(precip_annual),
             !is.na(temp_seasonality),
             !is.na(precip_seasonality),
             !is.na(Genus))
    
    # Merge variables by addition for annual model
    X1 <- scale(complete_data$temp_annual + complete_data$temp_seasonality)
    X2 <- scale(complete_data$precip_annual + complete_data$precip_seasonality)
  }
  
  # Common variables
  Y <- as.matrix(complete_data[[response_var]])
  X3 <- scale(complete_data$Sample_Age)
  X4 <- as.numeric(complete_data$Genus == "Oryza")
  
  cat("Sample size:", nrow(complete_data), "\n")
  
  # Perform 4-way variance partitioning
  vp <- varpart(Y, X1, X2, X3, X4)
  
  # --- Helpers ---------------------------------------------------------------
  # convert p-value to stars using your exact thresholds
  signif_stars <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("ns")
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    return("ns")
  }
  
  # safe RDA + ANOVA wrapper; returns list(F=..., p=...)
  test_rda <- function(formula_obj) {
    rda_obj <- tryCatch(rda(formula_obj), error = function(e) NULL)
    if (is.null(rda_obj)) return(list(F = NA, p = NA))
    an <- tryCatch(anova(rda_obj, permutations = 999), error = function(e) NULL)
    if (is.null(an)) return(list(F = NA, p = NA))
    return(list(F = an$F[1], p = an$`Pr(>F)`[1]))
  }
  # ---------------------------------------------------------------------------
  
  # --- Individual fractions --------------------------------------------------
  r_temp   <- test_rda(Y ~ X1 + Condition(X2 + X3 + X4))
  r_precip <- test_rda(Y ~ X2 + Condition(X1 + X3 + X4))
  r_age    <- test_rda(Y ~ X3 + Condition(X1 + X2 + X4))
  r_genus  <- test_rda(Y ~ X4 + Condition(X1 + X2 + X3))
  
  # --- Pairwise combinations (6) --------------------------------------------
  r_temp_precip  <- test_rda(Y ~ X1 + X2 + Condition(X3 + X4))
  r_temp_age     <- test_rda(Y ~ X1 + X3 + Condition(X2 + X4))
  r_temp_genus   <- test_rda(Y ~ X1 + X4 + Condition(X2 + X3))
  r_precip_age   <- test_rda(Y ~ X2 + X3 + Condition(X1 + X4))
  r_precip_genus <- test_rda(Y ~ X2 + X4 + Condition(X1 + X3))
  r_age_genus    <- test_rda(Y ~ X3 + X4 + Condition(X1 + X2))
  
  # --- Three-way combinations (4) -------------------------------------------
  r_temp_precip_age   <- test_rda(Y ~ X1 + X2 + X3 + Condition(X4))
  r_temp_precip_genus <- test_rda(Y ~ X1 + X2 + X4 + Condition(X3))
  r_temp_age_genus    <- test_rda(Y ~ X1 + X3 + X4 + Condition(X2))
  r_precip_age_genus  <- test_rda(Y ~ X2 + X3 + X4 + Condition(X1))
  
  # --- Full 4-way model -----------------------------------------------------
  r_full <- test_rda(Y ~ X1 + X2 + X3 + X4)
  
  # --- Print results to console ---------------------------------------------
  cat("\nSignificance testing of variance fractions using RDA:\n")
  cat(rep("-", 60), "\n")
  
  cat(sprintf("Individual fractions:\n"))
  cat(sprintf("Temp. [a]:   F = %s, p = %s %s\n",
              ifelse(is.na(r_temp$F), "NA", sprintf("%.3f", r_temp$F)),
              ifelse(is.na(r_temp$p), "NA", sprintf("%.3f", r_temp$p)),
              signif_stars(r_temp$p)))
  cat(sprintf("Precip. [b]: F = %s, p = %s %s\n",
              ifelse(is.na(r_precip$F), "NA", sprintf("%.3f", r_precip$F)),
              ifelse(is.na(r_precip$p), "NA", sprintf("%.3f", r_precip$p)),
              signif_stars(r_precip$p)))
  cat(sprintf("Age [c]:     F = %s, p = %s %s\n",
              ifelse(is.na(r_age$F), "NA", sprintf("%.3f", r_age$F)),
              ifelse(is.na(r_age$p), "NA", sprintf("%.3f", r_age$p)),
              signif_stars(r_age$p)))
  cat(sprintf("Genus [d]:   F = %s, p = %s %s\n",
              ifelse(is.na(r_genus$F), "NA", sprintf("%.3f", r_genus$F)),
              ifelse(is.na(r_genus$p), "NA", sprintf("%.3f", r_genus$p)),
              signif_stars(r_genus$p)))
  
  cat("\nPairwise combinations:\n")
  cat(sprintf("Temp+Precip:   F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_precip$F), "NA", sprintf("%.3f", r_temp_precip$F)),
              ifelse(is.na(r_temp_precip$p), "NA", sprintf("%.3f", r_temp_precip$p)),
              signif_stars(r_temp_precip$p)))
  cat(sprintf("Temp+Age:      F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_age$F), "NA", sprintf("%.3f", r_temp_age$F)),
              ifelse(is.na(r_temp_age$p), "NA", sprintf("%.3f", r_temp_age$p)),
              signif_stars(r_temp_age$p)))
  cat(sprintf("Temp+Genus:    F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_genus$F), "NA", sprintf("%.3f", r_temp_genus$F)),
              ifelse(is.na(r_temp_genus$p), "NA", sprintf("%.3f", r_temp_genus$p)),
              signif_stars(r_temp_genus$p)))
  cat(sprintf("Precip+Age:    F = %s, p = %s %s\n",
              ifelse(is.na(r_precip_age$F), "NA", sprintf("%.3f", r_precip_age$F)),
              ifelse(is.na(r_precip_age$p), "NA", sprintf("%.3f", r_precip_age$p)),
              signif_stars(r_precip_age$p)))
  cat(sprintf("Precip+Genus:  F = %s, p = %s %s\n",
              ifelse(is.na(r_precip_genus$F), "NA", sprintf("%.3f", r_precip_genus$F)),
              ifelse(is.na(r_precip_genus$p), "NA", sprintf("%.3f", r_precip_genus$p)),
              signif_stars(r_precip_genus$p)))
  cat(sprintf("Age+Genus:     F = %s, p = %s %s\n",
              ifelse(is.na(r_age_genus$F), "NA", sprintf("%.3f", r_age_genus$F)),
              ifelse(is.na(r_age_genus$p), "NA", sprintf("%.3f", r_age_genus$p)),
              signif_stars(r_age_genus$p)))
  
  cat("\nThree-way combinations:\n")
  cat(sprintf("Temp+Precip+Age:   F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_precip_age$F), "NA", sprintf("%.3f", r_temp_precip_age$F)),
              ifelse(is.na(r_temp_precip_age$p), "NA", sprintf("%.3f", r_temp_precip_age$p)),
              signif_stars(r_temp_precip_age$p)))
  cat(sprintf("Temp+Precip+Genus: F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_precip_genus$F), "NA", sprintf("%.3f", r_temp_precip_genus$F)),
              ifelse(is.na(r_temp_precip_genus$p), "NA", sprintf("%.3f", r_temp_precip_genus$p)),
              signif_stars(r_temp_precip_genus$p)))
  cat(sprintf("Temp+Age+Genus:    F = %s, p = %s %s\n",
              ifelse(is.na(r_temp_age_genus$F), "NA", sprintf("%.3f", r_temp_age_genus$F)),
              ifelse(is.na(r_temp_age_genus$p), "NA", sprintf("%.3f", r_temp_age_genus$p)),
              signif_stars(r_temp_age_genus$p)))
  cat(sprintf("Precip+Age+Genus:  F = %s, p = %s %s\n",
              ifelse(is.na(r_precip_age_genus$F), "NA", sprintf("%.3f", r_precip_age_genus$F)),
              ifelse(is.na(r_precip_age_genus$p), "NA", sprintf("%.3f", r_precip_age_genus$p)),
              signif_stars(r_precip_age_genus$p)))
  
  cat("\nOverall (4-way) model:\n")
  cat(sprintf("Full model:  F = %s, p = %s %s\n",
              ifelse(is.na(r_full$F), "NA", sprintf("%.3f", r_full$F)),
              ifelse(is.na(r_full$p), "NA", sprintf("%.3f", r_full$p)),
              signif_stars(r_full$p)))
  cat(rep("-", 60), "\n")
  
  # --- Plot Venn diagram ----------------------------------------------------
  plot(vp, Xnames = c("Temp.", "Precip.", "Age", "Genus"),
       bg = c("lightcoral", "lightblue", "lightgreen", "lightyellow"),
       digits = 3, cex = 0.8)
  
  # Add title
  title(main = paste("5' C>T Damage (Corrected) -", toupper(model_type), "Climate"))
  
  # Add Adj. R² + asterisks for full model significance
  full_model_lm <- lm(Y ~ X1 + X2 + X3 + X4)
  adj_r2 <- summary(full_model_lm)$adj.r.squared
  mtext(paste0("Adj. R² = ", sprintf("%.3f", adj_r2), " ", signif_stars(r_full$p)), 
        side = 3, line = 0.5, cex = 0.8)
  
  return(vp)
}

# Run Analysis for 5' C>T (Corrected)
cat("\n", rep("#", 80), "\n")
cat("VARIANCE PARTITIONING ANALYSIS: 5' C>T (Corrected)\n")
cat(rep("#", 80), "\n")
# Set up side-by-side plot layout
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
# Run collection climate model
vp_collection <- run_main_varpart_analysis(data_analysis, "X5P_DMG_POS1_Corrected", "collection")
# Run annual climate model
vp_annual <- run_main_varpart_analysis(data_analysis, "X5P_DMG_POS1_Corrected", "annual")
# Reset par
par(mfrow = c(1, 1))

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")

#------------------------------------------------------------------------------#
#  4.3 Regression analysis: 5' C>T (Corrected)  ~ Annual and Collection Temp   #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#               Combined Genera (single regression line)                       #
#------------------------------------------------------------------------------#

# Models for combined data
combined_annual_model <- lm(X5P_DMG_POS1_Corrected ~ temp_annual, data = data_analysis)
combined_coll_model <- lm(X5P_DMG_POS1_Corrected ~ Collection_Temp, data = data_analysis)

# Extract stats for combined models
combined_annual_r2 <- summary(combined_annual_model)$r.squared
combined_annual_p <- summary(combined_annual_model)$coefficients[2, 4]
combined_coll_r2 <- summary(combined_coll_model)$r.squared
combined_coll_p <- summary(combined_coll_model)$coefficients[2, 4]

# Plot 1: Annual temperature (combined regression)
Suppl_P6_5CT_Corrected_Temp <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_DMG_POS1_Corrected, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Annual Mean Temperature (°C)", y = "5' C>T Frequencies (Corrected)", tag = "a") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4,
           label = sprintf("R² = %.3f, p = %.3g", combined_annual_r2, combined_annual_p))

# Plot 2: Collection temperature (combined regression)
Suppl_P7_5CT_Corrected_Temp <- ggplot(data_analysis, aes(x = Collection_Temp, y = X5P_DMG_POS1_Corrected, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Temperature (°C)", y = "5' C>T Frequencies (Corrected)", tag = "b") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4,
           label = sprintf("R² = %.3f, p = %.3g", combined_coll_r2, combined_coll_p))

#------------------------------------------------------------------------------#
#               Separate Genera Analyses (genus-specific regression lines)    #
#------------------------------------------------------------------------------#

# Split data by genus for separate models
hordeum_data <- data_analysis %>% filter(Genus == "Hordeum")
oryza_data <- data_analysis %>% filter(Genus == "Oryza")

# Models for Hordeum
hordeum_annual_model <- lm(X5P_DMG_POS1_Corrected ~ temp_annual, data = hordeum_data)
hordeum_coll_model <- lm(X5P_DMG_POS1_Corrected ~ Collection_Temp, data = hordeum_data)

# Models for Oryza
oryza_annual_model <- lm(X5P_DMG_POS1_Corrected ~ temp_annual, data = oryza_data)
oryza_coll_model <- lm(X5P_DMG_POS1_Corrected ~ Collection_Temp, data = oryza_data)

# Extract stats for Hordeum
hordeum_annual_r2 <- summary(hordeum_annual_model)$r.squared
hordeum_annual_p <- summary(hordeum_annual_model)$coefficients[2, 4]
hordeum_coll_r2 <- summary(hordeum_coll_model)$r.squared
hordeum_coll_p <- summary(hordeum_coll_model)$coefficients[2, 4]

# Extract stats for Oryza
oryza_annual_r2 <- summary(oryza_annual_model)$r.squared
oryza_annual_p <- summary(oryza_annual_model)$coefficients[2, 4]
oryza_coll_r2 <- summary(oryza_coll_model)$r.squared
oryza_coll_p <- summary(oryza_coll_model)$coefficients[2, 4]

# Plot 3: Annual temperature (separate regressions by genus)
Suppl_P8_5CT_Corrected_Temp <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_DMG_POS1_Corrected, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Annual Mean Temperature (°C)", y = "5' C>T Frequencies (Corrected)", tag = "c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = c(1.5, 3), size = 3.5,
           label = c(sprintf("italic(Hordeum)~':  R²  =  %.3f,  p  =  %.3g'", hordeum_annual_r2, hordeum_annual_p),
                     sprintf("italic(Oryza)~':  R²  =  %.3f,  p  =  %.3g'", oryza_annual_r2, oryza_annual_p)),
           parse = TRUE)

# Plot 4: Collection temperature (separate regressions by genus)
Suppl_P9_5CT_Corrected_Temp <- ggplot(data_analysis, aes(x = Collection_Temp, y = X5P_DMG_POS1_Corrected, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Temperature (°C)", y = "5' C>T Frequencies (Corrected)", tag = "d") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = c(1.5, 3), size = 3.5,
           label = c(sprintf("italic(Hordeum)~':  R²  =  %.3f,  p  =  %.3g'", hordeum_coll_r2, hordeum_coll_p),
                     sprintf("italic(Oryza)~':  R²  =  %.3f,  p  =  %.3g'", oryza_coll_r2, oryza_coll_p)),
           parse = TRUE)
#------------------------------------------------------------------------------#
#               Create shared legend and combine plots                          #
#------------------------------------------------------------------------------#

# 1. Create a separate plot just for the legend
legend_plot <- Suppl_P6_5CT_Corrected_Temp + theme(legend.position = "bottom")
# 2. Extract just the legend
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
# 3. Remove legends from plots and add labels
Suppl_P6_5CT_Corrected_Temp <- Suppl_P6_5CT_Corrected_Temp + theme(legend.position = "none") + labs(tag = "a")
Suppl_P7_5CT_Corrected_Temp <- Suppl_P7_5CT_Corrected_Temp + theme(legend.position = "none") + labs(tag = "b")
Suppl_P8_5CT_Corrected_Temp <- Suppl_P8_5CT_Corrected_Temp + theme(legend.position = "none") + labs(tag = "c")
Suppl_P9_5CT_Corrected_Temp <- Suppl_P9_5CT_Corrected_Temp + theme(legend.position = "none") + labs(tag = "d")

# 4. Create the 2x2 grid with legend at bottom
grid.arrange(
  arrangeGrob(Suppl_P6_5CT_Corrected_Temp, Suppl_P7_5CT_Corrected_Temp, Suppl_P8_5CT_Corrected_Temp, Suppl_P9_5CT_Corrected_Temp, ncol = 2),
  legend,
  heights = c(10, 1),
  ncol = 1
)

#------------------------------------------------------------------------------#
#     Comparison of Raw vs Corrected 5' C>T Frequencies by Species            #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                    Statistical Testing by Species                            #
#------------------------------------------------------------------------------#
cat("\n=== Paired t-tests: Raw vs Corrected 5' C>T by Species ===\n\n")

# Perform paired t-test for each species
species_list <- unique(d_filtered$Species)
test_results <- data.frame()

for (sp in species_list) {
  sp_data <- d_filtered %>% filter(Species == sp)
  
  # Paired t-test (same samples, two conditions)
  test <- t.test(sp_data$X5P_DMG_POS1, sp_data$X5P_DMG_POS1_Corrected, 
                 paired = TRUE)
  
  # Calculate means and SDs
  mean_raw <- mean(sp_data$X5P_DMG_POS1, na.rm = TRUE)
  sd_raw <- sd(sp_data$X5P_DMG_POS1, na.rm = TRUE)
  mean_corrected <- mean(sp_data$X5P_DMG_POS1_Corrected, na.rm = TRUE)
  sd_corrected <- sd(sp_data$X5P_DMG_POS1_Corrected, na.rm = TRUE)
  mean_diff <- mean_raw - mean_corrected
  sd_diff <- sd(sp_data$X5P_DMG_POS1 - sp_data$X5P_DMG_POS1_Corrected, na.rm = TRUE)
  
  # Store results
  result_row <- data.frame(
    Species = sp,
    n = nrow(sp_data),
    Mean_Raw = mean_raw * 100,
    SD_Raw = sd_raw * 100,
    Mean_Corrected = mean_corrected * 100,
    SD_Corrected = sd_corrected * 100,
    Mean_Difference = mean_diff * 100,
    SD_Difference = sd_diff * 100,
    Percent_Reduction = (mean_diff / mean_raw) * 100,
    t_statistic = test$statistic,
    p_value = test$p.value,
    Significant = ifelse(test$p.value < 0.001, "***",
                         ifelse(test$p.value < 0.01, "**",
                                ifelse(test$p.value < 0.05, "*", "ns")))
  )
  
  test_results <- rbind(test_results, result_row)
  
  # Print results
  cat(sprintf("%-20s (n=%3d):\n", sp, nrow(sp_data)))
  cat(sprintf("  Raw:       %.3f ± %.3f%%\n", mean_raw*100, sd_raw*100))
  cat(sprintf("  Corrected: %.3f ± %.3f%%\n", mean_corrected*100, sd_corrected*100))
  cat(sprintf("  Difference: %.3f ± %.3f%% (%.1f%% reduction)\n", 
              mean_diff*100, sd_diff*100, (mean_diff/mean_raw)*100))
  cat(sprintf("  t = %.3f, p = %.2e %s\n\n", 
              test$statistic, test$p.value, 
              ifelse(test$p.value < 0.001, "***",
                     ifelse(test$p.value < 0.01, "**",
                            ifelse(test$p.value < 0.05, "*", "ns")))))
}

# Print summary table
cat("\n=== Summary Table ===\n")
print(test_results, row.names = FALSE)

# Add summary statistics as additional rows
cat("\n=== Summary Statistics ===\n")
summary_row_mean <- data.frame(
  Species = "MEAN",
  n = sum(test_results$n),
  Mean_Raw = mean(test_results$Mean_Raw),
  SD_Raw = NA,
  Mean_Corrected = mean(test_results$Mean_Corrected),
  SD_Corrected = NA,
  Mean_Difference = mean(test_results$Mean_Difference),
  SD_Difference = NA,
  Percent_Reduction = mean(test_results$Percent_Reduction),
  t_statistic = NA,
  p_value = NA,
  Significant = ""
)

summary_row_sd <- data.frame(
  Species = "SD",
  n = NA,
  Mean_Raw = sd(test_results$Mean_Raw),
  SD_Raw = NA,
  Mean_Corrected = sd(test_results$Mean_Corrected),
  SD_Corrected = NA,
  Mean_Difference = sd(test_results$Mean_Difference),
  SD_Difference = NA,
  Percent_Reduction = sd(test_results$Percent_Reduction),
  t_statistic = NA,
  p_value = NA,
  Significant = ""
)

# Combine with original results
test_results_with_summary <- rbind(test_results, summary_row_mean, summary_row_sd)

# Print with summary
print(test_results_with_summary, row.names = FALSE)

# Export single table
write.table(test_results_with_summary, 
            "5pCT_Raw_vs_Corrected_Statistical_Tests.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Table exported: 5pCT_Raw_vs_Corrected_Statistical_Tests.txt ===\n")
