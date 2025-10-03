#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#                     aDNA Damage through time - Script 06                     #
#                    Regression 5'C>T and Temperature models                   #
#------------------------------------------------------------------------------#

# Load required libraries
library(tidyverse)
library(terra)
library(ggplot2)
library(viridis)
library(gridExtra)
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
  mutate(Sample_Age = current_year - Collection_year)

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
#          Regression analysis: 5' C>T  ~ Annual and Collection Temp           #
#------------------------------------------------------------------------------#

# Define genus colors
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")

#------------------------------------------------------------------------------#
#               Combined Genera (single regression line)                       #
#------------------------------------------------------------------------------#

# Models for combined data
combined_annual_model <- lm(X5P_DMG_POS1 ~ temp_annual, data = data_analysis)
combined_coll_model <- lm(X5P_DMG_POS1 ~ Collection_Temp, data = data_analysis)

# Extract stats for combined models
combined_annual_r2 <- summary(combined_annual_model)$r.squared
combined_annual_p <- summary(combined_annual_model)$coefficients[2, 4]
combined_coll_r2 <- summary(combined_coll_model)$r.squared
combined_coll_p <- summary(combined_coll_model)$coefficients[2, 4]

# Plot 1: Annual temperature (combined regression)
P1 <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Annual Mean Temperature (°C)", y = "5' C>T Frequencies", tag = "a") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4,
           label = sprintf("R² = %.3f, p = %.3g", combined_annual_r2, combined_annual_p))

# Plot 2: Collection temperature (combined regression)
P2 <- ggplot(data_analysis, aes(x = Collection_Temp, y = X5P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Temperature (°C)", y = "5' C>T Frequencies", tag = "b") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4,
           label = sprintf("R² = %.3f, p = %.3g", combined_coll_r2, combined_coll_p))

#------------------------------------------------------------------------------#
#               Create shared legend and combine plots                          #
#------------------------------------------------------------------------------#

# Create a plot with legend on the bottom for extraction
legend_plot <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_DMG_POS1, color = Genus)) +
  geom_point() +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza"))),
                     name = "Genus") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm"))

# Extract legend using grid
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

shared_legend <- g_legend(legend_plot)

# Create the final plot using gridExtra
final_plot <- grid.arrange(
  arrangeGrob(P1, P2, ncol = 1, nrow = 2),
  shared_legend,
  ncol = 1,
  heights = c(10, 1)
)

# Save and display
grid.draw(final_plot)

#------------------------------------------------------------------------------#
#                             Print model summaries                            #
#------------------------------------------------------------------------------#

cat("=== SEPARATE GENUS MODELS ===\n")
cat("\nHordeum - Annual Temperature:\n")
print(summary(hordeum_annual_model))
cat("\nOryza - Annual Temperature:\n")
print(summary(oryza_annual_model))
cat("\nHordeum - Collection Temperature:\n")
print(summary(hordeum_coll_model))
cat("\nOryza - Collection Temperature:\n")
print(summary(oryza_coll_model))

cat("\n=== COMBINED MODELS ===\n")
cat("\nCombined - Annual Temperature:\n")
print(summary(combined_annual_model))
cat("\nCombined - Collection Temperature:\n")
print(summary(combined_coll_model))
