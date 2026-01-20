#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#                 aDNA Damage through time -  Script 07                        #
#                           SUPPLEMENTARY ANALYSES                             #
#                     Divergence by Species (Corrected)                        #
#------------------------------------------------------------------------------#

# Load libraries
library(dplyr)
library(ggplot2)

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


#------------------------------------------------------------------------------#
#                    Calculate divergence for each sample                      #
#------------------------------------------------------------------------------#
# For each sample
# 1. Misincorporation = Freq(C>T + G>A)
d_filtered$Deamination_rates <- as.numeric(d_filtered$X5P_DMG_POS1) + as.numeric(d_filtered$X3P_DMG_POS1)

# 2. Baseline sub rates = Freq(AllOthers)
d_filtered$Baseline_substitution_rates <- as.numeric(d_filtered$X5P_other_freq) + as.numeric(d_filtered$X3P_other_freq)

# 3. Corrected divergence = Freq(C>T + G>A) - Freq(AllOthers)
d_filtered$Corrected_Divergence <- d_filtered$Deamination_rates - d_filtered$Baseline_substitution_rates

# 4. Divergence effect on C>T = Freq(AllOthers) / Freq(C>T + G>A)
d_filtered$Divergence_Effect <- d_filtered$Baseline_substitution_rates / d_filtered$Deamination_rates
d_filtered$Divergence_Effect[is.infinite(d_filtered$Divergence_Effect)] <- NA
d_filtered$Divergence_Effect[is.nan(d_filtered$Divergence_Effect)] <- NA

#------------------------------------------------------------------------------#
#                    Summary by Species                                        #
#------------------------------------------------------------------------------#
species_summary <- d_filtered %>%
  group_by(Species) %>%
  summarise(
    n_samples = n(),
    Mean_Corrected_Divergence = mean(Corrected_Divergence, na.rm = TRUE),
    SD_Corrected_Divergence = sd(Corrected_Divergence, na.rm = TRUE),
    Mean_Divergence_Effect = mean(Divergence_Effect, na.rm = TRUE),
    SD_Divergence_Effect = sd(Divergence_Effect, na.rm = TRUE),
    .groups = "drop"
  )

print(species_summary, n = Inf)

#------------------------------------------------------------------------------#
#                    Boxplot                                                   #
#------------------------------------------------------------------------------#
# Define colors by genus
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")

P1 <- ggplot(d_filtered, aes(x = Species, y = Corrected_Divergence, fill = Genus)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = genus_colors,
                    labels = c(expression(italic("Hordeum")),
                              expression(italic("Oryza")))) +
  labs(x = NULL,
       y = "Corrected Divergence from Reference",
       title = "Reference genome divergence by species (corrected for deamination)",
       fill = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

print(P1)

#------------------------------------------------------------------------------#
#                    Export                                                    #
#------------------------------------------------------------------------------#
write.table(species_summary, 
            file = "Divergence_by_Species.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Table exported: Divergence_by_Species.txt ===\n")
