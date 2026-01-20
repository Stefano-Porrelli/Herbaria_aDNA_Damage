#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#            Patterns of aDNA Damage Through Time end Environments             #
#                      – lessons from herbarium specimens  -                   #
#                                                                              #
#                                   Script 05b                                 #
#                                                                              #
#                         SUPPLEMENTARY DATA ANALYSIS V-b:                     #
#                                                                              #
#                          CHELSA Climatic variables &                         #
#                     Variance partitioning (varpart/VEGAN)                    #
#                                                                              #
#                                                                              #
#                     - 3-way: Temperature, Precipitation, Age                 #
#                     - 4-way: Climate, Age, Genus, Herbarium                  #
#                      - Significance testing (rda + ANOVA)                    #
#------------------------------------------------------------------------------#

# Install required packages
# install.packages(c("tidyverse", "terra", "ggplot2", "viridis", "gridExtra", "vegan"))

# Load required libraries
library(tidyverse)
library(terra)
library(ggplot2)
library(viridis)
library(gridExtra)
library(vegan)

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

# Filter for Hordeum and Oryza genera only
d_genera <- d_filtered %>%
  filter(Genus %in% c("Hordeum", "Oryza"))

# Filter out rows with non-finite or zero values and calculate sample age
current_year <- 2025
d_clean <- d_genera %>%
  filter(!is.na(MEDIAN_SIZE) & MEDIAN_SIZE > 0,
         !is.na(Log_Mean),
         !is.na(X5P_DMG_POS1),
         !is.na(Lambda),
         !is.na(Collection_year),
         !is.na(Collection_month),
         !is.na(Herbarium)) %>%
  mutate(Sample_Age = current_year - Collection_year,
         # Create corrected 5' C>T value by subtracting non-deamination background
         X5P_DMG_POS1_Corrected = X5P_DMG_POS1 - X5P_other_freq,
         X3P_DMG_POS1_Corrected = X3P_DMG_POS1 - X3P_other_freq)

# Clean coordinates
clean_coords <- d_clean %>%
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

#------------------------------------------------------------------------------#
#                     Variance Partitioning analysis (VEGAN)                   #
#                               Response: Endogenous %                         #
#                     3-way: Temperature, Precipitation, Age                   # 
#------------------------------------------------------------------------------#

# --- Helpers ---------------------------------------------------------------
# convert p-value to stars
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

#------------------------------------------------------------------------------#
#                          COLLECTION Climate                                  #
#------------------------------------------------------------------------------#
Y      <- data_analysis$Endogenous_fraction
X1 <- data_analysis$Collection_Temp
X2 <- data_analysis$Collection_Precip
X3 <- data_analysis$Sample_Age

# 3-way variance partitioning (Temp + Precip + Age)
vp_collection <- varpart(Y,
                         ~ X1,
                         ~ X2,
                         ~ X3,
                         data = data_analysis)
print(vp_collection)

# --- Individual fractions
r_temp   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2 <- RsquareAdj(rda_full)$adj.r.squared

#------------------------------------------------------------------------------#
#                          ANNUAL Merged Climate                                #
#------------------------------------------------------------------------------#
X1 <- data_analysis$temp_annual + data_analysis$temp_seasonality
X2 <- data_analysis$precip_annual + data_analysis$precip_seasonality
X3 <- data_analysis$Sample_Age

vp_annual <- varpart(Y,
                     ~ X1,
                     ~ X2,
                     ~ X3,
                     data = data_analysis)
print(vp_annual)

# --- Individual fractions
r_temp_a   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip_a <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age_a    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip_a <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age_a    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age_a  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full_a <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full_a <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2_a <- RsquareAdj(rda_full_a)$adj.r.squared

#------------------------------------------------------------------------------#
#                          Side-by-side Venn diagrams                           #
#------------------------------------------------------------------------------#
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# --- Collection ---
plot(vp_collection, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2, 3), " ", signif_stars(r_full$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Endogenous % - COLLECTION Climate")

# --- Annual merged ---
plot(vp_annual, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2_a, 3), " ", signif_stars(r_full_a$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Endogenous % - ANNUAL Climate")

par(mfrow = c(1, 1))

#------------------------------------------------------------------------------#
#                       Print all fraction significances                        #
#------------------------------------------------------------------------------#
cat("\n--- COLLECTION Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp$F), "NA", sprintf("%.3f", r_temp$F)),
            ifelse(is.na(r_temp$p), "NA", sprintf("%.3f", r_temp$p)),
            signif_stars(r_temp$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip$F), "NA", sprintf("%.3f", r_precip$F)),
            ifelse(is.na(r_precip$p), "NA", sprintf("%.3f", r_precip$p)),
            signif_stars(r_precip$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age$F), "NA", sprintf("%.3f", r_age$F)),
            ifelse(is.na(r_age$p), "NA", sprintf("%.3f", r_age$p)),
            signif_stars(r_age$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip$F), "NA", sprintf("%.3f", r_temp_precip$F)),
            ifelse(is.na(r_temp_precip$p), "NA", sprintf("%.3f", r_temp_precip$p)),
            signif_stars(r_temp_precip$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age$F), "NA", sprintf("%.3f", r_temp_age$F)),
            ifelse(is.na(r_temp_age$p), "NA", sprintf("%.3f", r_temp_age$p)),
            signif_stars(r_temp_age$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age$F), "NA", sprintf("%.3f", r_precip_age$F)),
            ifelse(is.na(r_precip_age$p), "NA", sprintf("%.3f", r_precip_age$p)),
            signif_stars(r_precip_age$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full$F), "NA", sprintf("%.3f", r_full$F)),
            ifelse(is.na(r_full$p), "NA", sprintf("%.3f", r_full$p)),
            signif_stars(r_full$p)))
cat(rep("-", 60), "\n")

cat("\n--- ANNUAL Merged Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_a$F), "NA", sprintf("%.3f", r_temp_a$F)),
            ifelse(is.na(r_temp_a$p), "NA", sprintf("%.3f", r_temp_a$p)),
            signif_stars(r_temp_a$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_a$F), "NA", sprintf("%.3f", r_precip_a$F)),
            ifelse(is.na(r_precip_a$p), "NA", sprintf("%.3f", r_precip_a$p)),
            signif_stars(r_precip_a$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age_a$F), "NA", sprintf("%.3f", r_age_a$F)),
            ifelse(is.na(r_age_a$p), "NA", sprintf("%.3f", r_age_a$p)),
            signif_stars(r_age_a$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip_a$F), "NA", sprintf("%.3f", r_temp_precip_a$F)),
            ifelse(is.na(r_temp_precip_a$p), "NA", sprintf("%.3f", r_temp_precip_a$p)),
            signif_stars(r_temp_precip_a$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age_a$F), "NA", sprintf("%.3f", r_temp_age_a$F)),
            ifelse(is.na(r_temp_age_a$p), "NA", sprintf("%.3f", r_temp_age_a$p)),
            signif_stars(r_temp_age_a$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age_a$F), "NA", sprintf("%.3f", r_precip_age_a$F)),
            ifelse(is.na(r_precip_age_a$p), "NA", sprintf("%.3f", r_precip_age_a$p)),
            signif_stars(r_precip_age_a$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full_a$F), "NA", sprintf("%.3f", r_full_a$F)),
            ifelse(is.na(r_full_a$p), "NA", sprintf("%.3f", r_full_a$p)),
            signif_stars(r_full_a$p)))
cat(rep("-", 60), "\n")

#------------------------------------------------------------------------------#
#                     Variance Partitioning analysis (VEGAN)                   #
#                               Response: Fragment Size                        #
#                     3-way: Temperature, Precipitation, Age                   #
#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#
#                          COLLECTION Climate                                  #
#------------------------------------------------------------------------------#
Y      <- data_analysis$MEDIAN_SIZE
X1 <- data_analysis$Collection_Temp
X2 <- data_analysis$Collection_Precip
X3 <- data_analysis$Sample_Age

# 3-way variance partitioning (Temp + Precip + Age)
vp_collection <- varpart(Y,
                         ~ X1,
                         ~ X2,
                         ~ X3,
                         data = data_analysis)
print(vp_collection)

# --- Individual fractions
r_temp   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2 <- RsquareAdj(rda_full)$adj.r.squared

#------------------------------------------------------------------------------#
#                          ANNUAL Merged Climate                                #
#------------------------------------------------------------------------------#
X1 <- data_analysis$temp_annual + data_analysis$temp_seasonality
X2 <- data_analysis$precip_annual + data_analysis$precip_seasonality
X3 <- data_analysis$Sample_Age

vp_annual <- varpart(Y,
                     ~ X1,
                     ~ X2,
                     ~ X3,
                     data = data_analysis)
print(vp_annual)

# --- Individual fractions
r_temp_a   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip_a <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age_a    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip_a <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age_a    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age_a  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full_a <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full_a <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2_a <- RsquareAdj(rda_full_a)$adj.r.squared

#------------------------------------------------------------------------------#
#                          Side-by-side Venn diagrams                           #
#------------------------------------------------------------------------------#
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# --- Collection ---
plot(vp_collection, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2, 3), " ", signif_stars(r_full$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Fragment Size - COLLECTION Climate")

# --- Annual merged ---
plot(vp_annual, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2_a, 3), " ", signif_stars(r_full_a$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Fragment Size - ANNUAL Climate")

par(mfrow = c(1, 1))

#------------------------------------------------------------------------------#
#                       Print all fraction significances                        #
#------------------------------------------------------------------------------#
cat("\n--- COLLECTION Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp$F), "NA", sprintf("%.3f", r_temp$F)),
            ifelse(is.na(r_temp$p), "NA", sprintf("%.3f", r_temp$p)),
            signif_stars(r_temp$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip$F), "NA", sprintf("%.3f", r_precip$F)),
            ifelse(is.na(r_precip$p), "NA", sprintf("%.3f", r_precip$p)),
            signif_stars(r_precip$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age$F), "NA", sprintf("%.3f", r_age$F)),
            ifelse(is.na(r_age$p), "NA", sprintf("%.3f", r_age$p)),
            signif_stars(r_age$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip$F), "NA", sprintf("%.3f", r_temp_precip$F)),
            ifelse(is.na(r_temp_precip$p), "NA", sprintf("%.3f", r_temp_precip$p)),
            signif_stars(r_temp_precip$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age$F), "NA", sprintf("%.3f", r_temp_age$F)),
            ifelse(is.na(r_temp_age$p), "NA", sprintf("%.3f", r_temp_age$p)),
            signif_stars(r_temp_age$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age$F), "NA", sprintf("%.3f", r_precip_age$F)),
            ifelse(is.na(r_precip_age$p), "NA", sprintf("%.3f", r_precip_age$p)),
            signif_stars(r_precip_age$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full$F), "NA", sprintf("%.3f", r_full$F)),
            ifelse(is.na(r_full$p), "NA", sprintf("%.3f", r_full$p)),
            signif_stars(r_full$p)))
cat(rep("-", 60), "\n")

cat("\n--- ANNUAL Merged Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_a$F), "NA", sprintf("%.3f", r_temp_a$F)),
            ifelse(is.na(r_temp_a$p), "NA", sprintf("%.3f", r_temp_a$p)),
            signif_stars(r_temp_a$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_a$F), "NA", sprintf("%.3f", r_precip_a$F)),
            ifelse(is.na(r_precip_a$p), "NA", sprintf("%.3f", r_precip_a$p)),
            signif_stars(r_precip_a$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age_a$F), "NA", sprintf("%.3f", r_age_a$F)),
            ifelse(is.na(r_age_a$p), "NA", sprintf("%.3f", r_age_a$p)),
            signif_stars(r_age_a$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip_a$F), "NA", sprintf("%.3f", r_temp_precip_a$F)),
            ifelse(is.na(r_temp_precip_a$p), "NA", sprintf("%.3f", r_temp_precip_a$p)),
            signif_stars(r_temp_precip_a$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age_a$F), "NA", sprintf("%.3f", r_temp_age_a$F)),
            ifelse(is.na(r_temp_age_a$p), "NA", sprintf("%.3f", r_temp_age_a$p)),
            signif_stars(r_temp_age_a$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age_a$F), "NA", sprintf("%.3f", r_precip_age_a$F)),
            ifelse(is.na(r_precip_age_a$p), "NA", sprintf("%.3f", r_precip_age_a$p)),
            signif_stars(r_precip_age_a$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full_a$F), "NA", sprintf("%.3f", r_full_a$F)),
            ifelse(is.na(r_full_a$p), "NA", sprintf("%.3f", r_full_a$p)),
            signif_stars(r_full_a$p)))
cat(rep("-", 60), "\n")

#------------------------------------------------------------------------------#
#                     Variance Partitioning analysis (VEGAN)                   #
#                       Response: 5' C>T damage                                #
#                     3-way: Temperature, Precipitation, Age                   # 
#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#
#                          COLLECTION Climate                                  #
#------------------------------------------------------------------------------#
Y      <- data_analysis$X5P_DMG_POS1_Corrected
X1 <- data_analysis$Collection_Temp
X2 <- data_analysis$Collection_Precip
X3 <- data_analysis$Sample_Age

# 3-way variance partitioning (Temp + Precip + Age)
vp_collection <- varpart(Y,
                         ~ X1,
                         ~ X2,
                         ~ X3,
                         data = data_analysis)
print(vp_collection)

# --- Individual fractions
r_temp   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2 <- RsquareAdj(rda_full)$adj.r.squared

#------------------------------------------------------------------------------#
#                          ANNUAL Merged Climate                                #
#------------------------------------------------------------------------------#
X1 <- data_analysis$temp_annual + data_analysis$temp_seasonality
X2 <- data_analysis$precip_annual + data_analysis$precip_seasonality
X3 <- data_analysis$Sample_Age

vp_annual <- varpart(Y,
                     ~ X1,
                     ~ X2,
                     ~ X3,
                     data = data_analysis)
print(vp_annual)

# --- Individual fractions
r_temp_a   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip_a <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age_a    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip_a <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age_a    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age_a  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full_a <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full_a <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2_a <- RsquareAdj(rda_full_a)$adj.r.squared

#------------------------------------------------------------------------------#
#                          Side-by-side Venn diagrams                           #
#------------------------------------------------------------------------------#
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# --- Collection ---
plot(vp_collection, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2, 3), " ", signif_stars(r_full$p)),
      side = 3, line = 0.1, cex = 0.9)
title("5' C>T Damage - COLLECTION Climate")

# --- Annual merged ---
plot(vp_annual, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2_a, 3), " ", signif_stars(r_full_a$p)),
      side = 3, line = 0.1, cex = 0.9)
title("5' C>T Damage - ANNUAL Climate")

par(mfrow = c(1, 1))

#------------------------------------------------------------------------------#
#                       Print all fraction significances                        #
#------------------------------------------------------------------------------#
cat("\n--- COLLECTION Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp$F), "NA", sprintf("%.3f", r_temp$F)),
            ifelse(is.na(r_temp$p), "NA", sprintf("%.3f", r_temp$p)),
            signif_stars(r_temp$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip$F), "NA", sprintf("%.3f", r_precip$F)),
            ifelse(is.na(r_precip$p), "NA", sprintf("%.3f", r_precip$p)),
            signif_stars(r_precip$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age$F), "NA", sprintf("%.3f", r_age$F)),
            ifelse(is.na(r_age$p), "NA", sprintf("%.3f", r_age$p)),
            signif_stars(r_age$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip$F), "NA", sprintf("%.3f", r_temp_precip$F)),
            ifelse(is.na(r_temp_precip$p), "NA", sprintf("%.3f", r_temp_precip$p)),
            signif_stars(r_temp_precip$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age$F), "NA", sprintf("%.3f", r_temp_age$F)),
            ifelse(is.na(r_temp_age$p), "NA", sprintf("%.3f", r_temp_age$p)),
            signif_stars(r_temp_age$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age$F), "NA", sprintf("%.3f", r_precip_age$F)),
            ifelse(is.na(r_precip_age$p), "NA", sprintf("%.3f", r_precip_age$p)),
            signif_stars(r_precip_age$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full$F), "NA", sprintf("%.3f", r_full$F)),
            ifelse(is.na(r_full$p), "NA", sprintf("%.3f", r_full$p)),
            signif_stars(r_full$p)))
cat(rep("-", 60), "\n")

cat("\n--- ANNUAL Merged Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_a$F), "NA", sprintf("%.3f", r_temp_a$F)),
            ifelse(is.na(r_temp_a$p), "NA", sprintf("%.3f", r_temp_a$p)),
            signif_stars(r_temp_a$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_a$F), "NA", sprintf("%.3f", r_precip_a$F)),
            ifelse(is.na(r_precip_a$p), "NA", sprintf("%.3f", r_precip_a$p)),
            signif_stars(r_precip_a$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age_a$F), "NA", sprintf("%.3f", r_age_a$F)),
            ifelse(is.na(r_age_a$p), "NA", sprintf("%.3f", r_age_a$p)),
            signif_stars(r_age_a$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip_a$F), "NA", sprintf("%.3f", r_temp_precip_a$F)),
            ifelse(is.na(r_temp_precip_a$p), "NA", sprintf("%.3f", r_temp_precip_a$p)),
            signif_stars(r_temp_precip_a$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age_a$F), "NA", sprintf("%.3f", r_temp_age_a$F)),
            ifelse(is.na(r_temp_age_a$p), "NA", sprintf("%.3f", r_temp_age_a$p)),
            signif_stars(r_temp_age_a$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age_a$F), "NA", sprintf("%.3f", r_precip_age_a$F)),
            ifelse(is.na(r_precip_age_a$p), "NA", sprintf("%.3f", r_precip_age_a$p)),
            signif_stars(r_precip_age_a$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full_a$F), "NA", sprintf("%.3f", r_full_a$F)),
            ifelse(is.na(r_full_a$p), "NA", sprintf("%.3f", r_full_a$p)),
            signif_stars(r_full_a$p)))
cat(rep("-", 60), "\n")

#------------------------------------------------------------------------------#
#                     Variance Partitioning analysis (VEGAN)                   #
#                               Response: Lambda                               #
#                     3-way: Temperature, Precipitation, Age                   # 
#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#
#                          COLLECTION Climate                                  #
#------------------------------------------------------------------------------#
Y      <- data_analysis$Lambda
X1 <- data_analysis$Collection_Temp
X2 <- data_analysis$Collection_Precip
X3 <- data_analysis$Sample_Age

# 3-way variance partitioning (Temp + Precip + Age)
vp_collection <- varpart(Y,
                         ~ X1,
                         ~ X2,
                         ~ X3,
                         data = data_analysis)
print(vp_collection)

# --- Individual fractions
r_temp   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2 <- RsquareAdj(rda_full)$adj.r.squared

#------------------------------------------------------------------------------#
#                          ANNUAL Merged Climate                                #
#------------------------------------------------------------------------------#
X1 <- data_analysis$temp_annual + data_analysis$temp_seasonality
X2 <- data_analysis$precip_annual + data_analysis$precip_seasonality
X3 <- data_analysis$Sample_Age

vp_annual <- varpart(Y,
                     ~ X1,
                     ~ X2,
                     ~ X3,
                     data = data_analysis)
print(vp_annual)

# --- Individual fractions
r_temp_a   <- test_rda(Y ~ X1 + Condition(X2 + X3))
r_precip_a <- test_rda(Y ~ X2 + Condition(X1 + X3))
r_age_a    <- test_rda(Y ~ X3 + Condition(X1 + X2))

# --- Pairwise combinations
r_temp_precip_a <- test_rda(Y ~ X1 + X2 + Condition(X3))
r_temp_age_a    <- test_rda(Y ~ X1 + X3 + Condition(X2))
r_precip_age_a  <- test_rda(Y ~ X2 + X3 + Condition(X1))

# --- Full 3-way model
r_full_a <- test_rda(Y ~ X1 + X2 + X3)

# --- RDA for full model Adj. R²
rda_full_a <- rda(Y ~ X1 + X2 + X3, data = data_analysis)
adjR2_a <- RsquareAdj(rda_full_a)$adj.r.squared

#------------------------------------------------------------------------------#
#                          Side-by-side Venn diagrams                           #
#------------------------------------------------------------------------------#
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# --- Collection ---
plot(vp_collection, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2, 3), " ", signif_stars(r_full$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Lambda - COLLECTION Climate")

# --- Annual merged ---
plot(vp_annual, digits = 2,
     bg = c("tomato", "skyblue", "palegreen3"),
     Xnames = c("Temp", "Precip", "Age"))
mtext(paste0("Adj R² = ", round(adjR2_a, 3), " ", signif_stars(r_full_a$p)),
      side = 3, line = 0.1, cex = 0.9)
title("Lambda - ANNUAL Climate")

par(mfrow = c(1, 1))

#------------------------------------------------------------------------------#
#                       Print all fraction significances                        #
#------------------------------------------------------------------------------#
cat("\n--- COLLECTION Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp$F), "NA", sprintf("%.3f", r_temp$F)),
            ifelse(is.na(r_temp$p), "NA", sprintf("%.3f", r_temp$p)),
            signif_stars(r_temp$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip$F), "NA", sprintf("%.3f", r_precip$F)),
            ifelse(is.na(r_precip$p), "NA", sprintf("%.3f", r_precip$p)),
            signif_stars(r_precip$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age$F), "NA", sprintf("%.3f", r_age$F)),
            ifelse(is.na(r_age$p), "NA", sprintf("%.3f", r_age$p)),
            signif_stars(r_age$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip$F), "NA", sprintf("%.3f", r_temp_precip$F)),
            ifelse(is.na(r_temp_precip$p), "NA", sprintf("%.3f", r_temp_precip$p)),
            signif_stars(r_temp_precip$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age$F), "NA", sprintf("%.3f", r_temp_age$F)),
            ifelse(is.na(r_temp_age$p), "NA", sprintf("%.3f", r_temp_age$p)),
            signif_stars(r_temp_age$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age$F), "NA", sprintf("%.3f", r_precip_age$F)),
            ifelse(is.na(r_precip_age$p), "NA", sprintf("%.3f", r_precip_age$p)),
            signif_stars(r_precip_age$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full$F), "NA", sprintf("%.3f", r_full$F)),
            ifelse(is.na(r_full$p), "NA", sprintf("%.3f", r_full$p)),
            signif_stars(r_full$p)))
cat(rep("-", 60), "\n")

cat("\n--- ANNUAL Merged Climate ---\n")
cat(sprintf("Temp [a]:   F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_a$F), "NA", sprintf("%.3f", r_temp_a$F)),
            ifelse(is.na(r_temp_a$p), "NA", sprintf("%.3f", r_temp_a$p)),
            signif_stars(r_temp_a$p)))
cat(sprintf("Precip [b]: F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_a$F), "NA", sprintf("%.3f", r_precip_a$F)),
            ifelse(is.na(r_precip_a$p), "NA", sprintf("%.3f", r_precip_a$p)),
            signif_stars(r_precip_a$p)))
cat(sprintf("Age [c]:    F = %s, p = %s %s\n",
            ifelse(is.na(r_age_a$F), "NA", sprintf("%.3f", r_age_a$F)),
            ifelse(is.na(r_age_a$p), "NA", sprintf("%.3f", r_age_a$p)),
            signif_stars(r_age_a$p)))
cat("\nPairwise combinations:\n")
cat(sprintf("Temp+Precip: F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_precip_a$F), "NA", sprintf("%.3f", r_temp_precip_a$F)),
            ifelse(is.na(r_temp_precip_a$p), "NA", sprintf("%.3f", r_temp_precip_a$p)),
            signif_stars(r_temp_precip_a$p)))
cat(sprintf("Temp+Age:    F = %s, p = %s %s\n",
            ifelse(is.na(r_temp_age_a$F), "NA", sprintf("%.3f", r_temp_age_a$F)),
            ifelse(is.na(r_temp_age_a$p), "NA", sprintf("%.3f", r_temp_age_a$p)),
            signif_stars(r_temp_age_a$p)))
cat(sprintf("Precip+Age:  F = %s, p = %s %s\n",
            ifelse(is.na(r_precip_age_a$F), "NA", sprintf("%.3f", r_precip_age_a$F)),
            ifelse(is.na(r_precip_age_a$p), "NA", sprintf("%.3f", r_precip_age_a$p)),
            signif_stars(r_precip_age_a$p)))
cat("\nOverall (3-way) model:\n")
cat(sprintf("Full model:  F = %s, p = %s %s\n",
            ifelse(is.na(r_full_a$F), "NA", sprintf("%.3f", r_full_a$F)),
            ifelse(is.na(r_full_a$p), "NA", sprintf("%.3f", r_full_a$p)),
            signif_stars(r_full_a$p)))
cat(rep("-", 60), "\n")





#------------------------------------------------------------------------------#
#                     Variance Partitioning analysis (VEGAN)                   #
#             Responses: Endogenous %, Fragment size, 5'C>T, Lambda            #
#                     4-way: Climate, Age, Genus, Herbarium                    # 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                  Supplementary Variance Partitioning Function                #
#------------------------------------------------------------------------------#
run_supplementary_varpart_analysis <- function(data, response_var, model_type) {
  
  cat("\n", rep("=", 70), "\n")
  cat("SUPPLEMENTARY ANALYSIS:", response_var, "-", toupper(model_type), "CLIMATE + HERBARIUM\n")
  cat(rep("=", 70), "\n")
  
  # Filter data and prepare variables based on model type
  if (model_type == "collection") {
    complete_data <- data %>%
      filter(!is.na(!!sym(response_var)),
             !is.na(Sample_Age),
             !is.na(Collection_Temp),
             !is.na(Collection_Precip),
             !is.na(Genus),
             !is.na(Herbarium))
    
    # Merged climate for collection model
    X1 <- scale(complete_data$Collection_Temp + complete_data$Collection_Precip)
    
  } else if (model_type == "annual") {
    complete_data <- data %>%
      filter(!is.na(!!sym(response_var)),
             !is.na(Sample_Age),
             !is.na(temp_annual),
             !is.na(precip_annual),
             !is.na(temp_seasonality),
             !is.na(precip_seasonality),
             !is.na(Genus),
             !is.na(Herbarium))
    
    # Merged climate for annual model (all 4 variables)
    X1 <- scale(complete_data$temp_annual + complete_data$precip_annual + 
               complete_data$temp_seasonality + complete_data$precip_seasonality)
  }
  
  # Common variables
  Y <- as.matrix(complete_data[[response_var]])
  X2 <- scale(complete_data$Sample_Age)
  X3 <- as.numeric(complete_data$Genus == "Oryza")
  X4 <- as.numeric(as.factor(complete_data$Herbarium))
  
  cat("Sample size:", nrow(complete_data), "\n")
  cat("Number of herbaria:", length(unique(complete_data$Herbarium)), "\n")
  
  # 4-way variance partitioning
  vp <- varpart(Y, X1, X2, X3, X4)
  
  # --- Helpers ---------------------------------------------------------------
  signif_stars <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("ns")
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    return("ns")
  }
  
  test_rda <- function(formula_obj) {
    rda_obj <- tryCatch(rda(formula_obj), error = function(e) NULL)
    if (is.null(rda_obj)) return(list(F = NA, p = NA))
    an <- tryCatch(anova(rda_obj, permutations = 999), error = function(e) NULL)
    if (is.null(an)) return(list(F = NA, p = NA))
    return(list(F = an$F[1], p = an$`Pr(>F)`[1]))
  }
  # ---------------------------------------------------------------------------
  
  # --- Individual fractions --------------------------------------------------
  r_climate   <- test_rda(Y ~ X1 + Condition(X2 + X3 + X4))
  r_age       <- test_rda(Y ~ X2 + Condition(X1 + X3 + X4))
  r_genus     <- test_rda(Y ~ X3 + Condition(X1 + X2 + X4))
  r_herbarium <- test_rda(Y ~ X4 + Condition(X1 + X2 + X3))
  
  # --- Pairwise combinations (6) --------------------------------------------
  r_climate_age      <- test_rda(Y ~ X1 + X2 + Condition(X3 + X4))
  r_climate_genus    <- test_rda(Y ~ X1 + X3 + Condition(X2 + X4))
  r_climate_herbarium<- test_rda(Y ~ X1 + X4 + Condition(X2 + X3))
  r_age_genus        <- test_rda(Y ~ X2 + X3 + Condition(X1 + X4))
  r_age_herbarium    <- test_rda(Y ~ X2 + X4 + Condition(X1 + X3))
  r_genus_herbarium  <- test_rda(Y ~ X3 + X4 + Condition(X1 + X2))
  
  # --- Three-way combinations (4) -------------------------------------------
  r_climate_age_genus       <- test_rda(Y ~ X1 + X2 + X3 + Condition(X4))
  r_climate_age_herbarium   <- test_rda(Y ~ X1 + X2 + X4 + Condition(X3))
  r_climate_genus_herbarium <- test_rda(Y ~ X1 + X3 + X4 + Condition(X2))
  r_age_genus_herbarium     <- test_rda(Y ~ X2 + X3 + X4 + Condition(X1))
  
  # --- Full 4-way model -----------------------------------------------------
  r_full <- test_rda(Y ~ X1 + X2 + X3 + X4)
  
  # --- Print results --------------------------------------------------------
  cat("\nSignificance testing of variance fractions using RDA:\n")
  cat(rep("-", 60), "\n")
  
  cat("Individual fractions:\n")
  cat(sprintf("Climate [a]:   F = %s, p = %s %s\n", r_climate$F, r_climate$p, signif_stars(r_climate$p)))
  cat(sprintf("Age [b]:       F = %s, p = %s %s\n", r_age$F, r_age$p, signif_stars(r_age$p)))
  cat(sprintf("Genus [c]:     F = %s, p = %s %s\n", r_genus$F, r_genus$p, signif_stars(r_genus$p)))
  cat(sprintf("Herbarium [d]: F = %s, p = %s %s\n", r_herbarium$F, r_herbarium$p, signif_stars(r_herbarium$p)))
  
  cat("\nPairwise combinations:\n")
  cat(sprintf("Climate+Age:          F = %s, p = %s %s\n", r_climate_age$F, r_climate_age$p, signif_stars(r_climate_age$p)))
  cat(sprintf("Climate+Genus:        F = %s, p = %s %s\n", r_climate_genus$F, r_climate_genus$p, signif_stars(r_climate_genus$p)))
  cat(sprintf("Climate+Herbarium:    F = %s, p = %s %s\n", r_climate_herbarium$F, r_climate_herbarium$p, signif_stars(r_climate_herbarium$p)))
  cat(sprintf("Age+Genus:            F = %s, p = %s %s\n", r_age_genus$F, r_age_genus$p, signif_stars(r_age_genus$p)))
  cat(sprintf("Age+Herbarium:        F = %s, p = %s %s\n", r_age_herbarium$F, r_age_herbarium$p, signif_stars(r_age_herbarium$p)))
  cat(sprintf("Genus+Herbarium:      F = %s, p = %s %s\n", r_genus_herbarium$F, r_genus_herbarium$p, signif_stars(r_genus_herbarium$p)))
  
  cat("\nThree-way combinations:\n")
  cat(sprintf("Climate+Age+Genus:          F = %s, p = %s %s\n", r_climate_age_genus$F, r_climate_age_genus$p, signif_stars(r_climate_age_genus$p)))
  cat(sprintf("Climate+Age+Herbarium:      F = %s, p = %s %s\n", r_climate_age_herbarium$F, r_climate_age_herbarium$p, signif_stars(r_climate_age_herbarium$p)))
  cat(sprintf("Climate+Genus+Herbarium:    F = %s, p = %s %s\n", r_climate_genus_herbarium$F, r_climate_genus_herbarium$p, signif_stars(r_climate_genus_herbarium$p)))
  cat(sprintf("Age+Genus+Herbarium:        F = %s, p = %s %s\n", r_age_genus_herbarium$F, r_age_genus_herbarium$p, signif_stars(r_age_genus_herbarium$p)))
  
  cat("\nOverall (4-way) model:\n")
  cat(sprintf("Full model:    F = %s, p = %s %s\n", r_full$F, r_full$p, signif_stars(r_full$p)))
  cat(rep("-", 60), "\n")
  
  # --- Plot Venn diagram ----------------------------------------------------
  plot(vp, Xnames = c("Climate", "Age", "Genus", "Herbarium"),
       bg = c("lightcoral", "lightgreen", "lightyellow", "lightgray"),
       digits = 3, cex = 0.8)
  
  response_name <- switch(response_var,
                          "X5P_DMG_POS1" = "5' C>T Damage",
                          "MEDIAN_SIZE" = "Fragment Size", 
                          "Lambda" = "Lambda",
                          "Endogenous_fraction" = "Endogenous Fraction",
                          response_var)
  title(main = paste(response_name, "-", toupper(model_type)))
  
  full_model_lm <- lm(Y ~ X1 + X2 + X3 + X4)
  adj_r2 <- summary(full_model_lm)$adj.r.squared
  mtext(paste0("Adj. R² = ", sprintf("%.3f", adj_r2), " ", signif_stars(r_full$p)),
        side = 3, line = 0.5, cex = 0.8)
  
  return(vp)
}

#------------------------------------------------------------------------------#
#                        Run Supplementary Analyses                            #
#------------------------------------------------------------------------------#

response_vars <- c("X5P_DMG_POS1", "MEDIAN_SIZE", "Lambda", "Endogenous_fraction")
response_names <- c("5' C>T Damage", "Fragment Size", "Lambda", "Endogenous Fraction")

for (i in seq_along(response_vars)) {
  
  cat("\n", rep("#", 80), "\n")
  cat("SUPPLEMENTARY ANALYSIS FOR:", response_names[i], "\n")
  cat(rep("#", 80), "\n")
  
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  
  vp_supp_collection <- run_supplementary_varpart_analysis(data_analysis, response_vars[i], "collection")
  vp_supp_annual <- run_supplementary_varpart_analysis(data_analysis, response_vars[i], "annual")

  par(mfrow = c(1, 1))
}
