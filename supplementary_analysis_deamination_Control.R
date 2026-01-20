#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#  SUPPLEMENTARY MISINCORPORATIONS ANALYSIS                                    #
#  Extract all misincorporation types from MapDamage output                    #
#  Calculate "other" average                                                   #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# 1: Extract misincorporations from MapDamage files
#------------------------------------------------------------------------------#

# Function to extract misincorporations from misincorporation.txt
extract_other_misincorp <- function(file_path) {
  # Read the file
  data <- read.table(file_path, header = TRUE, comment.char = "#", 
                     stringsAsFactors = FALSE)
  
  # Filter for position 1 at both 5' and 3' ends
  pos1_5p <- data %>% filter(End == "5p", Pos == 1)
  pos1_3p <- data %>% filter(End == "3p", Pos == 1)
  
  # Calculate for 5' end (ALL misincorporations including C>T for comparison)
  if (nrow(pos1_5p) > 0) {
    totals_5p <- pos1_5p %>%
      summarise(
        Total_A = sum(A),
        Total_T = sum(T),
        Total_G = sum(G),
        Total_C = sum(C),
        # ALL misincorporations at 5' end
        # NOTE: R converts > to . in column names when reading tables
        C_to_T_5p = sum(C.T),
        G_to_A_5p = sum(G.A),
        A_to_G_5p = sum(A.G),
        A_to_C_5p = sum(A.C),
        A_to_T_5p = sum(A.T),
        T_to_C_5p = sum(T.C),
        T_to_G_5p = sum(T.G),
        T_to_A_5p = sum(T.A),
        G_to_C_5p = sum(G.C),
        G_to_T_5p = sum(G.T),
        C_to_G_5p = sum(C.G),
        C_to_A_5p = sum(C.A)
      )
  } else {
    return(NULL)
  }
  
  # Calculate for 3' end (ALL misincorporations including G>A for comparison)
  if (nrow(pos1_3p) > 0) {
    totals_3p <- pos1_3p %>%
      summarise(
        # ALL misincorporations at 3' end
        # NOTE: R converts > to . in column names when reading tables
        G_to_A_3p = sum(G.A),
        C_to_T_3p = sum(C.T),
        A_to_G_3p = sum(A.G),
        A_to_C_3p = sum(A.C),
        A_to_T_3p = sum(A.T),
        T_to_C_3p = sum(T.C),
        T_to_G_3p = sum(T.G),
        T_to_A_3p = sum(T.A),
        G_to_C_3p = sum(G.C),
        G_to_T_3p = sum(G.T),
        C_to_G_3p = sum(C.G),
        C_to_A_3p = sum(C.A)
      )
  } else {
    totals_3p <- data.frame()
  }
  
  # Combine 5' and 3' results
  result <- cbind(totals_5p, totals_3p)
  
  # Calculate frequencies (misincorp count / total source nucleotide)
  result <- result %>%
    mutate(
      # 5' end frequencies - DEAMINATION (calculated from misincorporation.txt)
      C_to_T_5p_calc_freq = C_to_T_5p / Total_C,
      G_to_A_5p_calc_freq = G_to_A_5p / Total_G,
      # 5' end frequencies - NON-DEAMINATION
      A_to_G_5p_freq = A_to_G_5p / Total_A,
      A_to_C_5p_freq = A_to_C_5p / Total_A,
      A_to_T_5p_freq = A_to_T_5p / Total_A,
      T_to_C_5p_freq = T_to_C_5p / Total_T,
      T_to_G_5p_freq = T_to_G_5p / Total_T,
      T_to_A_5p_freq = T_to_A_5p / Total_T,
      G_to_C_5p_freq = G_to_C_5p / Total_G,
      G_to_T_5p_freq = G_to_T_5p / Total_G,
      C_to_G_5p_freq = C_to_G_5p / Total_C,
      C_to_A_5p_freq = C_to_A_5p / Total_C,
      # 3' end frequencies - DEAMINATION (calculated from misincorporation.txt)
      G_to_A_3p_calc_freq = G_to_A_3p / Total_G,
      C_to_T_3p_calc_freq = C_to_T_3p / Total_C,
      # 3' end frequencies - NON-DEAMINATION
      A_to_G_3p_freq = A_to_G_3p / Total_A,
      A_to_C_3p_freq = A_to_C_3p / Total_A,
      A_to_T_3p_freq = A_to_T_3p / Total_A,
      T_to_C_3p_freq = T_to_C_3p / Total_T,
      T_to_G_3p_freq = T_to_G_3p / Total_T,
      T_to_A_3p_freq = T_to_A_3p / Total_T,
      G_to_C_3p_freq = G_to_C_3p / Total_G,
      G_to_T_3p_freq = G_to_T_3p / Total_G,
      C_to_G_3p_freq = C_to_G_3p / Total_C,
      C_to_A_3p_freq = C_to_A_3p / Total_C
    )
  
  return(result)
}

#------------------------------------------------------------------------------#
# Find all sample directories
#------------------------------------------------------------------------------#
sample_dirs <- unique(dirname(list.files(pattern = "misincorporation\\.txt$", 
                                         recursive = TRUE, 
                                         full.names = TRUE)))

#------------------------------------------------------------------------------#
# Process all samples
#------------------------------------------------------------------------------#
all_data <- data.frame()

for (sample_dir in sample_dirs) {
  sample_id <- basename(sample_dir)
  
  # File paths
  misinc_file <- file.path(sample_dir, "misincorporation.txt")
  ct_file <- file.path(sample_dir, "5pCtoT_freq.txt")
  ga_file <- file.path(sample_dir, "3pGtoA_freq.txt")
  
  # Check all files exist
  if (!file.exists(misinc_file) | !file.exists(ct_file) | !file.exists(ga_file)) {
    cat("  Skipping", sample_id, "- missing files\n")
    next
  }
  
  # Extract 5' C>T from dedicated file (position 1)
  # NOTE: R converts > to . in column names
  ct_data <- read.table(ct_file, header = TRUE)
  ct_freq <- ct_data[ct_data$pos == 1, 2]  # Second column (X5pC.T)
  
  # Extract 3' G>A from dedicated file (position 1)
  ga_data <- read.table(ga_file, header = TRUE)
  ga_freq <- ga_data[ga_data$pos == 1, 2]  # Second column (X3pG.A)
  
  # Extract all other misincorporations
  other_misinc <- extract_other_misincorp(misinc_file)
  
  if (!is.null(other_misinc)) {
    # Combine everything
    sample_data <- data.frame(
      Sample = sample_id,
      # From pre-calculated files (MapDamage standard output)
      C_to_T_5p_precalc = ct_freq,
      G_to_A_3p_precalc = ga_freq
    )
    # Add calculated values from misincorporation.txt
    sample_data <- cbind(sample_data, other_misinc)
    
    all_data <- rbind(all_data, sample_data)
  }
}

cat("\nProcessed", nrow(all_data), "samples successfully\n")

#------------------------------------------------------------------------------#
# PART 2: Calculate "other" substitutions averages
#------------------------------------------------------------------------------#

# Calculate "other" for 5' end (average of all non-deamination types)
all_data <- all_data %>%
  mutate(
    other_5p_freq = (A_to_G_5p_freq + A_to_C_5p_freq + A_to_T_5p_freq +
                     T_to_C_5p_freq + T_to_G_5p_freq + T_to_A_5p_freq +
                     G_to_C_5p_freq + G_to_T_5p_freq +
                     C_to_G_5p_freq + C_to_A_5p_freq) / 10,
    
    # Calculate "other" for 3' end (average of all non-deamination types)
    other_3p_freq = (A_to_G_3p_freq + A_to_C_3p_freq + A_to_T_3p_freq +
                     T_to_C_3p_freq + T_to_G_3p_freq + T_to_A_3p_freq +
                     G_to_C_3p_freq + G_to_T_3p_freq +
                     C_to_G_3p_freq + C_to_A_3p_freq) / 10
  )

# Save final results
write.csv(all_data, "all_misincorporations_with_other.csv", row.names = FALSE)

#------------------------------------------------------------------------------#
#  Merge misincorporation data with main metadata                             #
#  Match: Sample (misincorp) = Sample (metadata)                              #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Read data files
#------------------------------------------------------------------------------#
cat("Reading data files...\n")

# Read misincorporation data
misinc_data <- read.csv("all_misincorporations_with_other.csv", 
                        stringsAsFactors = FALSE)
cat("  Misincorporation data: ", nrow(misinc_data), "rows\n")

# Read main metadata
metadata <- read.delim("aDNA_damage_screening_MAIN.txt")
metadata <- as.data.frame(metadata)
cat("  Metadata: ", nrow(metadata), "rows\n")

#------------------------------------------------------------------------------#
# Check column names
#------------------------------------------------------------------------------#
cat("\nChecking key columns...\n")
cat("  Misincorporation key column: 'Sample'\n")
cat("  Metadata key column: 'Sample'\n")

if (!"Sample" %in% colnames(misinc_data)) {
  stop("ERROR: 'Sample' column not found in misincorporation data!")
}

if (!"Sample" %in% colnames(metadata)) {
  stop("ERROR: 'Sample' column not found in metadata!")
}

#------------------------------------------------------------------------------#
# Check for matching IDs
#------------------------------------------------------------------------------#
cat("\nChecking ID matches...\n")

# Get unique IDs from each dataset
misinc_ids <- unique(misinc_data$Sample)
metadata_ids <- unique(metadata$Sample)

cat("  Unique samples in misincorporation data:", length(misinc_ids), "\n")
cat("  Unique samples in metadata:", length(metadata_ids), "\n")

# Check overlap
matches <- sum(misinc_ids %in% metadata_ids)
cat("  Samples present in BOTH datasets:", matches, "\n")

# Samples in misincorp but NOT in metadata
missing_in_metadata <- misinc_ids[!misinc_ids %in% metadata_ids]
if (length(missing_in_metadata) > 0) {
  cat("\n  WARNING:", length(missing_in_metadata), 
      "samples in misincorp data but NOT in metadata:\n")
  cat("    ", head(missing_in_metadata, 10), "\n")
  if (length(missing_in_metadata) > 10) {
    cat("    ... and", length(missing_in_metadata) - 10, "more\n")
  }
}

# Samples in metadata but NOT in misincorp
missing_in_misinc <- metadata_ids[!metadata_ids %in% misinc_ids]
if (length(missing_in_misinc) > 0) {
  cat("\n  WARNING:", length(missing_in_misinc), 
      "samples in metadata but NOT in misincorp data:\n")
  cat("    ", head(missing_in_misinc, 10), "\n")
  if (length(missing_in_misinc) > 10) {
    cat("    ... and", length(missing_in_misinc) - 10, "more\n")
  }
}

#------------------------------------------------------------------------------#
# Merge datasets - select only 3 columns from misincorp data
#------------------------------------------------------------------------------#
cat("\n", rep("=", 80), "\n", sep = "")
cat("PERFORMING MERGE\n")
cat(rep("=", 80), "\n", sep = "")

# Select only the 3 columns we want from misincorporation data
misinc_subset <- misinc_data %>%
  select(Sample, 
         G_to_A_3p_precalc, 
         other_5p_freq, 
         other_3p_freq)

# Merge with metadata (keeps all metadata columns + adds 3 new ones)
merged_data <- metadata %>%
  left_join(misinc_subset, by = "Sample")

# Rename the new columns
merged_data <- merged_data %>%
  rename(
    `3P_DMG_POS1` = G_to_A_3p_precalc,
    `5P_other_freq` = other_5p_freq,
    `3P_other_freq` = other_3p_freq
  )

cat("\nMerge complete!\n")
cat("  Merged dataset: ", nrow(merged_data), "rows\n")
cat("  Columns: ", ncol(merged_data), "\n")

# Check how many rows have misincorporation data
rows_with_misinc <- sum(!is.na(merged_data$`3P_DMG_POS1`))
cat("  Rows with misincorporation data:", rows_with_misinc, "\n")
cat("  Rows without misincorporation data:", nrow(merged_data) - rows_with_misinc, "\n")

#------------------------------------------------------------------------------#
# Save merged data
#------------------------------------------------------------------------------#
cat("\n", rep("=", 80), "\n", sep = "")
cat("SAVING RESULTS\n")
cat(rep("=", 80), "\n", sep = "")

# Save as tab-delimited file (all original columns + 3 new ones)
write.table(merged_data, "aDNA_damage_screening_MAIN.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)
cat("Saved: aDNA_damage_screening_MAIN.txt\n")


#==========================================================================================#


#------------------------------------------------------------------------------#
#            Supplementary Analyses: non-Deamination mismatch controls         #
#                                      Age                                     #
#------------------------------------------------------------------------------#
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
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
#------------------------------------------------------------------------------#
#               1 - Regression analysis: 5' C>T  ~ Sample Age                  #
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
P1_5p_DMG_pos1_Genus <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_DMG_POS1, color = Genus)) +
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
P1_5p_DMG_pos1_Genus

#------------------------------------------------------------------------------#
#               2 - Regression analysis: 3' C>T  ~ Sample Age                  #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 3pCT ~ sample age - by Genera
hordeum_model <- lm(X3P_DMG_POS1 ~ Sample_Age,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X3P_DMG_POS1 ~ Sample_Age,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
P2_3p_DMG_pos1_Genus <- ggplot(d_filtered, aes(x = Sample_Age, y = X3P_DMG_POS1, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01), 
                      limits = c(0, 0.055)) +  # Adjust y-axis with proper increments
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "3' G>A Frequencies",
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
P2_3p_DMG_pos1_Genus


#------------------------------------------------------------------------------#
#               3 - Regression analysis: 5' C>T  ~ 3' G>A                      #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ 3pGA - by Genera
hordeum_model <- lm(X5P_DMG_POS1 ~ X3P_DMG_POS1,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X5P_DMG_POS1 ~ X3P_DMG_POS1,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
P3_5p_vs_3p_Genus <- ggplot(d_filtered, aes(x = X3P_DMG_POS1, y = X5P_DMG_POS1, color = Genus)) +
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
  annotate("text", size = 3.5, x = 0.035, y = 0.0035,
           label = sprintf("italic(Hordeum): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           hordeum_r2, hordeum_p),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 0.035, y = 0.000005,
           label = sprintf("italic(Oryza): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           oryza_r2, oryza_p),
           parse = TRUE, hjust = 0)
# print plot
P3_5p_vs_3p_Genus

#------------------------------------------------------------------------------#
#               4 - Regression analysis: 5' Others  ~ Sample Age               #
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
P4_5p_Others_vs_age <- ggplot(d_filtered, aes(x = Sample_Age, y = X5P_other_freq, color = Genus)) +
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
P4_5p_Others_vs_age

#------------------------------------------------------------------------------#
#               5 - Regression analysis: 3' Others  ~ Sample Age                  #
#                                Divided by Genera                             #
#------------------------------------------------------------------------------#
# Regression: 5pCT ~ sample age - by Genera
hordeum_model <- lm(X3P_other_freq ~ Sample_Age,
                    data = filter(d_filtered, Genus == "Hordeum"))
oryza_model <- lm(X3P_other_freq ~ Sample_Age,
                  data = filter(d_filtered, Genus == "Oryza"))
# Extract R-squared and p-values for each model
hordeum_r2 <- summary(hordeum_model)$r.squared
hordeum_p <- summary(hordeum_model)$coefficients[2, 4]
oryza_r2 <- summary(oryza_model)$r.squared
oryza_p <- summary(oryza_model)$coefficients[2, 4]
P5_3p_Others_vs_age <- ggplot(d_filtered, aes(x = Sample_Age, y = X3P_other_freq, color = Genus)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Genus)) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Sample Age (years)",
       y = "3' Other Frequencies",
       title = NULL,
       color = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
    legend.position = "bottom") +
  # Add R-squared and p-value annotations for both groups separately
  annotate("text", size = 3.5, x = 130, y = 0.0012,
           label = sprintf("italic(Hordeum): 
             italic(R)^2 == %.5f~~italic(p) == %.3g",
                           hordeum_r2, hordeum_p),
           parse = TRUE, hjust = 0) +
  annotate("text", size = 3.5, x = 130, y = 0.000005,
           label = sprintf("italic(Oryza): 
             italic(R)^2 == %.3f~~italic(p) == %.3g",
                           oryza_r2, oryza_p),
           parse = TRUE, hjust = 0)
# print plot
P5_3p_Others_vs_age


#------------------------------------------------------------------------------#
# Print main plots together with labels                                       #
#------------------------------------------------------------------------------#
# 1. Create a separate plot just for the legend
legend_plot <- P1_5p_DMG_pos1_Genus + theme(legend.position = "bottom")
# 2. Extract just the legend
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
# 3. Remove legends from plots and add labels
P1_5p_DMG_pos1_Genus <- P1_5p_DMG_pos1_Genus + theme(legend.position = "none") + labs(tag = "a")
P2_3p_DMG_pos1_Genus <- P2_3p_DMG_pos1_Genus + theme(legend.position = "none") + labs(tag = "b")
P4_5p_Others_vs_age <- P4_5p_Others_vs_age + theme(legend.position = "none") + labs(tag = "c")
P5_3p_Others_vs_age <- P5_3p_Others_vs_age + theme(legend.position = "none") + labs(tag = "d")

# 4. Create the 2x2 grid 
grid.arrange(
  arrangeGrob(P1_5p_DMG_pos1_Genus, P2_3p_DMG_pos1_Genus, P4_5p_Others_vs_age, P5_3p_Others_vs_age, ncol = 2),
  legend,
  heights = c(10, 1),
  ncol = 1
)

#------------------------------------------------------------------------------#
#            Supplementary Analyses: non-Deamination mismatch controls         #
#                                      Age                                     #
#------------------------------------------------------------------------------#
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
#          Regression analysis: 5' Others  ~ Annual and Collection Temp        #
#------------------------------------------------------------------------------#

# Define genus colors
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")

#------------------------------------------------------------------------------#
#               Combined Genera (single regression line)                       #
#------------------------------------------------------------------------------#

# Models for combined data
combined_annual_model <- lm(X5P_other_freq ~ temp_annual, data = data_analysis)
combined_coll_model <- lm(X5P_other_freq ~ Collection_Temp, data = data_analysis)

# Extract stats for combined models
combined_annual_r2 <- summary(combined_annual_model)$r.squared
combined_annual_p <- summary(combined_annual_model)$coefficients[2, 4]
combined_coll_r2 <- summary(combined_coll_model)$r.squared
combined_coll_p <- summary(combined_coll_model)$coefficients[2, 4]

# Plot 1: Annual temperature (combined regression)
P1 <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_other_freq, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Annual Mean Temperature (°C)", y = "5' Other Frequencies", tag = "a") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 4,
           label = sprintf("R² = %.3f, p = %.3g", combined_annual_r2, combined_annual_p))

# Plot 2: Collection temperature (combined regression)
P2 <- ggplot(data_analysis, aes(x = Collection_Temp, y = X5P_other_freq, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Temperature (°C)", y = "5' Other Frequencies", tag = "b") +
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
hordeum_annual_model <- lm(X5P_other_freq ~ temp_annual, data = hordeum_data)
hordeum_coll_model <- lm(X5P_other_freq ~ Collection_Temp, data = hordeum_data)

# Models for Oryza
oryza_annual_model <- lm(X5P_other_freq ~ temp_annual, data = oryza_data)
oryza_coll_model <- lm(X5P_other_freq ~ Collection_Temp, data = oryza_data)

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
P3 <- ggplot(data_analysis, aes(x = temp_annual, y = X5P_other_freq, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Annual Mean Temperature (°C)", y = "5' Other Frequencies", tag = "c") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = c(1.5, 3), size = 3.5,
           label = c(sprintf("italic(Hordeum)~':  R²  =  %.3f,  p  =  %.3g'", hordeum_annual_r2, hordeum_annual_p),
                    sprintf("italic(Oryza)~':  R²  =  %.3f,  p  =  %.3g'", oryza_annual_r2, oryza_annual_p)),
           parse = TRUE)

# Plot 4: Collection temperature (separate regressions by genus)
P4 <- ggplot(data_analysis, aes(x = Collection_Temp, y = X5P_other_freq, color = Genus)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = genus_colors,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  labs(x = "Collection Temperature (°C)", y = "5' Other Frequencies", tag = "d") +
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
legend_plot <- P1 + theme(legend.position = "bottom")
# 2. Extract just the legend
legend <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
# 3. Remove legends from plots and add labels
P1 <- P1 + theme(legend.position = "none") + labs(tag = "a")
P2 <- P2 + theme(legend.position = "none") + labs(tag = "b")
P3 <- P3 + theme(legend.position = "none") + labs(tag = "c")
P4 <- P4 + theme(legend.position = "none") + labs(tag = "d")

# 4. Create the 2x2 grid with legend at bottom
grid.arrange(
  arrangeGrob(P1, P2, P3, P4, ncol = 2),
  legend,
  heights = c(10, 1),
  ncol = 1
)
