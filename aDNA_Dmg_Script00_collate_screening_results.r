#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#            Patterns of aDNA Damage Through Time end Environments             #
#                      – lessons from herbarium specimens  -                   #
#                                                                              #
#                                   Script 00                                  #
#                                                                              #
#                               DATA PREPARATION:                              #
#                                                                              #
#   Merges species-specific screening directories, then collates results       #
#   from aDNA screening generated with Latorre et.al, 2020                     #
#     "Plant aDNA pipeline" (https://doi.org/10.1002/cppb.20121) protocol:     #
#               https://gitlab.com/smlatorreo/plant-adna-pipeline              #
#                                                                              #
#                        - Merge screening directories                         #
#                        - Collate screening results                           #
#   - Calculate lambda from exponential decay of fragment length distribution  #
#               - Plot exponential decay of fragment length                    #
#           - Calculate log-mean from fragment length distribution             #
#                                                                              #                                                                              #
# Note: sample_metadata.txt can be provided in a tab-separated file, ensure    #
#       names under the column "Sample" match the name of the samples          #
#       used during screening.                                                 #
#------------------------------------------------------------------------------#

# Set warning handling and messages
options(warn = 1)  # Print warnings as they occur
print("Starting script execution...")

# Install required packages, if needed uncomment below
# install.packages(c("dplyr", "readr", "MASS", "ggplot2"))

# Load necessary libraries
library(readr)
library(dplyr)
library(MASS)
library(ggplot2)

#------------------------------------------------------------------------------#
#                        STEP 1: Merge Screening Directories                  #
#------------------------------------------------------------------------------#

# Define screening directories to merge
screening_base_dir <- "./Screening_results"
species_dirs <- c(
  "Hordeum_spontaneum_RBGK",
  "Hordeum_vulgare_RBGK",
  "Oryza_Americas_UoN_KAUST",
  "Oryza_grandiglumis_SIHUS",
  "Oryza_latifolia_NYBG",
  "Oryza_rufipogon_RBGK"
)

subdirs_to_merge <- c(
  "2_trimmed_merged",
  "3_quality_control",
  "4_mapping",
  "5_aDNA_characteristics",
  "6_AMBER"
)

output_merged_dir <- "./Screening_results_collated"

# =================== Function to create directory structure ==================#
create_directory_structure <- function(output_dir, subdirs) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  for (subdir in subdirs) {
    subdir_path <- file.path(output_dir, subdir)
    if (!dir.exists(subdir_path)) {
      dir.create(subdir_path, recursive = TRUE)
    }
  }
}

# =================== Function to copy directory contents =====================#
copy_directory_contents <- function(source_dir, dest_dir) {
  files <- list.files(source_dir, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) return(0)
  
  files_copied <- 0
  for (file in files) {
    filename <- basename(file)
    dest_file <- file.path(dest_dir, filename)
    
    if (file.info(file)$isdir) {
      if (!dir.exists(dest_file)) dir.create(dest_file, recursive = TRUE)
      sub_files_copied <- copy_directory_contents(file, dest_file)
      files_copied <- files_copied + sub_files_copied
    } else {
      file.copy(file, dest_file, overwrite = FALSE)
      files_copied <- files_copied + 1
    }
  }
  return(files_copied)
}

# =================== Merge screening directories =============================#
message("\n=== Merging screening directories ===")
create_directory_structure(output_merged_dir, subdirs_to_merge)

total_files_copied <- 0
for (species in species_dirs) {
  message("Processing: ", species)
  species_path <- file.path(screening_base_dir, species)
  
  if (!dir.exists(species_path)) {
    warning("Species directory not found: ", species_path)
    next
  }
  
  for (subdir in subdirs_to_merge) {
    source_path <- file.path(species_path, subdir)
    dest_path <- file.path(output_merged_dir, subdir)
    
    if (dir.exists(source_path)) {
      files_copied <- copy_directory_contents(source_path, dest_path)
      total_files_copied <- total_files_copied + files_copied
    }
  }
}

message("Merge complete! Total files copied: ", total_files_copied)
message("======================================\n")

#------------------------------------------------------------------------------#
#                    STEP 2: Collate Screening Results                        #
#------------------------------------------------------------------------------#

# Define the base directories (pointing to merged directory)
base_dir <- output_merged_dir
mapdamage_dir <- file.path(base_dir, "5_aDNA_characteristics")
mapping_dir <- file.path(base_dir, "4_mapping")

#------------------------------------------------------------------------------#
#                                   Helper Functions                           #
#------------------------------------------------------------------------------#

# =================== Function to process flagstat files ======================#
process_flagstat <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(c(NA, NA))
  }
  lines <- readLines(file_path)
  mapped_line <- lines[grep("mapped \\(", lines)[1]]
  if (length(mapped_line) == 0) {
    warning("No mapped reads line found in flagstat file")
    return(c(NA, NA))
  }
  # Extract number of mapped reads
  no_merged <- as.numeric(strsplit(mapped_line, " ")[[1]][1])
  # Extract endogenous percentage
  percentage_match <- regexpr("\\([0-9.]+%", mapped_line)
  if (percentage_match > 0) {
    endogenous <- as.numeric(gsub("[^0-9.]", "",
                                  substr(mapped_line, percentage_match + 1,
                                         percentage_match +
                                           attr(percentage_match,
                                                "match.length") - 1)))
  } else {
    warning("No percentage found in mapped reads line")
    endogenous <- NA
  }
  return(c(endogenous, no_merged))
}

# =================== Function to extract duplication rate ======================#
get_duplication_rate <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NA)
  }
  lines <- readLines(file_path)
  dup_line <- tail(lines, 1)
  dup_rate <- as.numeric(gsub("Duplication Rate: ", "", dup_line))
  return(dup_rate)
}

# ========= Function to calculate lambda (damage fraction per site)===============#
calculate_lambda <- function(lengths, frequencies) {
  # Max fragment length to consider (exclude longer fragments)
  max_length <- 110
  # Create data frame for all data points
  all_data <- data.frame(
    length = lengths,
    frequency = frequencies
  ) %>%
    filter(frequency > 0)  # Remove zeros before log transformation
  # Find the peak using all data points
  peak_idx <- which.max(all_data$frequency)
  
  if (length(peak_idx) == 0 || peak_idx == nrow(all_data)) {
    warning("No clear peak found in distribution")
    peak_length <- min(all_data$length)
  } else {
    peak_length <- all_data$length[peak_idx]
  }
  # Create regression data - Only use points after peak and up to max_length
  regression_data <- all_data %>%
    filter(length >= peak_length) %>%  # Only after peak
    filter(length <= max_length) %>%   # Only up to max_length
    mutate(log_freq = log(frequency))
  # Check if we have enough points for regression
  if (nrow(regression_data) < 3) {
    warning("Not enough points for lambda calculation between peak and max_length")
    return(list(lambda = NA, r_squared = NA, p_value = NA))
  }
  # Perform linear regression
  model <- try(lm(log_freq ~ length, data = regression_data), silent = TRUE)
  if (inherits(model, "try-error") || is.null(model)) {
    warning("Failed to fit linear model for lambda calculation")
    return(list(lambda = NA, r_squared = NA, p_value = NA))
  }
  # Lambda is the negative of the slope
  lambda <- -coef(model)[2]
  r_squared <- summary(model)$r.squared
  # Calculate p-value for lambda (slope coefficient)
  p_value <- summary(model)$coefficients[2, 4]
  return(list(
    lambda = lambda,
    r_squared = r_squared,
    p_value = p_value
  ))
}

# ==================== Function to visualize the lambda plots ===========================#
plot_lambda_fit <- function(lengths, frequencies, lambda_results, sample_name, max_length = 110) {
  # Create data frame (show all but only use the ones that show exponential decay for fitting)
  all_data <- data.frame(
    length = lengths,
    frequency = frequencies
  ) %>% 
    filter(frequency > 0) %>%
    mutate(log_freq = log(frequency))
  
  # Find the peak for highlighting (using all data)
  peak_idx <- which.max(all_data$frequency)
  peak_length <- all_data$length[peak_idx]
  
  # Categorize data points for visualization - simplify to just "Used for fit" vs "Excluded"
  all_data$category <- "Excluded"
  all_data$category[all_data$length >= peak_length & all_data$length <= max_length] <- "Used for fit"
  
  # Create the regression line data
  if (!is.na(lambda_results$lambda)) {
    # Extract slope and intercept from the model
    slope <- -lambda_results$lambda
    
    # Find the exact data used for fitting (points after peak and <= max_length)
    fit_data <- all_data[all_data$category == "Used for fit",]
    
    # We need to recalculate the intercept for plotting
    intercept <- mean(fit_data$log_freq) - slope * mean(fit_data$length)
    
    # Generate prediction line
    x_range <- seq(min(fit_data$length), max(fit_data$length), length.out = 110)
    pred_line <- data.frame(
      length = x_range,
      log_freq = intercept + slope * x_range
    )
  } else {
    pred_line <- NULL
  }
  
  # Create plot with all data points
  p <- ggplot(all_data, aes(x = length, y = log_freq)) +
    geom_point(aes(color = category), alpha = 0.7) +
    scale_color_manual(values = c("grey60", "forestgreen"),
                       name = "Data Points") +
    # Add marker lines
    geom_vline(xintercept = peak_length, linetype = "dashed", color = "red") +
    geom_vline(xintercept = max_length, linetype = "dashed", color = "blue") +
    # Add labels
    annotate("text", x = peak_length + 5, y = max(all_data$log_freq), 
             label = "Peak", color = "red", hjust = 0) +
    annotate("text", x = max_length - 5, y = max(all_data$log_freq) - 0.5, 
             label = "110bp", color = "blue", hjust = 1) +
    labs(
      title = paste("Exponential fragment length distribution:", sample_name),
      subtitle = sprintf("Lambda = %.4f, R² = %.4f, p = %.6f", 
                         lambda_results$lambda, lambda_results$r_squared, lambda_results$p_value),
      x = "Fragment Length",
      y = "Log(Frequency)"
    ) +
    theme_minimal()
  
  # Add regression line if available
  if (!is.null(pred_line)) {
    p <- p + geom_line(data = pred_line, 
                       aes(x = length, y = log_freq), 
                       color = "red", size = 1)
  }
  
  return(p)
}

# ================ Function to analyze C-to-T damage decay pattern ======================#
analyze_dmg_decay <- function(file_path, max_pos = 20) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(list(
      decay_rate = NA,
      N_value = NA,
      t_value = NA,
      p_value = NA,
      df = NA,
      r_squared = NA
    ))
  }
  
  # Read damage data
  data <- try({
    read.table(file_path, header = TRUE, check.names = FALSE)
  }, silent = TRUE)
  
  if (inherits(data, "try-error") || nrow(data) < max_pos || !"5pC>T" %in% names(data)) {
    warning(paste("Invalid data format or insufficient data in:", file_path))
    return(list(
      decay_rate = NA,
      N_value = NA,
      t_value = NA,
      p_value = NA,
      df = NA,
      r_squared = NA
    ))
  }
  
  # Extract position and frequency data for first max_pos positions
  pos_data <- data.frame(
    position = 1:max_pos,
    freq = data[1:max_pos, "5pC>T"]
  )
  
  # Remove any NA or zero values that would cause problems with nls
  pos_data <- pos_data %>% 
    filter(!is.na(freq)) %>%
    filter(freq > 0)
  
  # Check if we have enough data points
  if (nrow(pos_data) < 5) {
    warning("Not enough valid data points for exponential fitting")
    return(list(
      decay_rate = NA,
      N_value = NA,
      t_value = NA,
      p_value = NA,
      df = NA,
      r_squared = NA
    ))
  }
  
  # Initial parameter estimates (starting values for nls)
  init_N <- max(pos_data$freq)
  init_rate <- 0.1  # Initial guess for decay rate
  
  # Fit exponential decay model using nls: y ~ N*exp(-rate*x)
  model <- try(
    nls(freq ~ N * exp(-rate * position), 
        data = pos_data,
        start = list(N = init_N, rate = init_rate),
        control = nls.control(maxiter = 100)),
    silent = TRUE
  )
  
  # Check if model fitting was successful
  if (inherits(model, "try-error") || is.null(model)) {
    warning("Failed to fit exponential model to damage data")
    return(list(
      decay_rate = NA,
      N_value = NA,
      t_value = NA,
      p_value = NA,
      df = NA,
      r_squared = NA
    ))
  }
  
  # Extract model parameters
  params <- summary(model)
  N_value <- params$coefficients["N", "Estimate"]
  decay_rate <- params$coefficients["rate", "Estimate"]
  
  # Get t-value and degrees of freedom for rate parameter
  t_value <- params$coefficients["rate", "t value"]
  df <- params$df[2]
  
  # Calculate p-value using one-sided t-test (we expect rate to be positive)
  p_value <- pt(t_value, df, lower.tail = FALSE)
  
  # Calculate R-squared (1 - RSS/TSS)
  fitted_values <- predict(model)
  residuals <- pos_data$freq - fitted_values
  RSS <- sum(residuals^2)
  TSS <- sum((pos_data$freq - mean(pos_data$freq))^2)
  r_squared <- 1 - RSS/TSS
  
  return(list(
    decay_rate = decay_rate,
    N_value = N_value,
    t_value = t_value,
    p_value = p_value,
    df = df,
    r_squared = r_squared,
    data = pos_data,  # Return the data for plotting
    model = model     # Return the model for prediction
  ))
}

# ================ Function to plot damage decay analysis ==============================#
plot_dmg_decay <- function(decay_results, file_name = "Unnamed Sample") {
  if (is.na(decay_results$decay_rate) || is.null(decay_results$data) || is.null(decay_results$model)) {
    warning("Cannot create plot with invalid model or missing data")
    return(NULL)
  }
  
  # Create data frame for plotting
  plot_data <- decay_results$data
  
  # Generate prediction line
  x_smooth <- seq(min(plot_data$position), max(plot_data$position), length.out = 100)
  y_pred <- decay_results$N_value * exp(-decay_results$decay_rate * x_smooth)
  pred_line <- data.frame(position = x_smooth, freq = y_pred)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = position, y = freq)) +
    geom_point(color = "blue", size = 3, alpha = 0.7) +
    geom_line(data = pred_line, aes(x = position, y = freq), 
              color = "red", size = 1) +
    labs(
      title = paste("C>T Damage Decay Pattern:", file_name),
      subtitle = sprintf(
        "R² = %.4f, p = %.6f\nModel: freq ~ %.4f * exp(-%.4f * position)",
        decay_results$r_squared,
        decay_results$p_value,
        decay_results$N_value,
        decay_results$decay_rate
      ),
      x = "Position from 5' End",
      y = "C>T Substitution Frequency"
    ) +
    theme_minimal() +
    ylim(0, max(plot_data$freq) * 1.1) # Add a little space above points
  
  return(p)
}

#================= Function to process length distribution and calculate metrics ==============#
process_lgdistribution <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(list(
      log_stats = c(NA, NA),
      lambda_stats = c(NA, NA, NA),
      median_size = NA,
      raw_distribution = NULL
    ))
  }
  
  # Read the length distribution data
  lgdist_data <- read_delim(file_path, delim = "\t", 
                            comment = "#", col_names = TRUE) %>%
    setNames(trimws(names(.)))
  
  # Create expanded lengths vector for calculations
  expanded_lengths <- rep(lgdist_data$Length, lgdist_data$Occurences)
  
  # Fit lognormal distribution
  fit <- try(fitdistr(expanded_lengths, "lognormal"), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    warning(paste("Failed to fit lognormal distribution for:", file_path))
    log_mean <- NA
    log_median <- NA
  } else {
    log_mean <- fit$estimate['meanlog']
    log_median <- exp(fit$estimate['meanlog'])
  }
  
  # Calculate median size
  median_size <- median(expanded_lengths)
  
  # Calculate lambda with 110bp limit
  lambda_results <- calculate_lambda(lgdist_data$Length, lgdist_data$Occurences)
  
  # Generate plot
  sample_name <- basename(dirname(file_path))
  plot <- try(plot_lambda_fit(lgdist_data$Length, lgdist_data$Occurences, 
                              lambda_results, sample_name, max_length = 110), silent = TRUE)
  
  if (!inherits(plot, "try-error") && !is.null(plot)) {
    plot_file <- file.path(dirname(dirname(file_path)), 
                           paste0(sample_name, "_lambda_fit.pdf"))
    try(ggsave(plot_file, plot, width = 8, height = 6), silent = TRUE)
  }
  
  return(list(
    log_stats = c(log_mean, log_median),
    lambda_stats = c(lambda_results$lambda, lambda_results$r_squared, lambda_results$p_value),
    median_size = median_size,
    raw_distribution = lgdist_data
  ))
}

#============== Function to extract damage metrics and analyze damage decay =============#
extract_damage_metrics <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(list(
      positions_1_3 = c(NA, NA, NA),
      decay_metrics = list(
        decay_rate = NA,
        N_value = NA,
        t_value = NA,
        p_value = NA,
        r_squared = NA
      )
    ))
  }
  
  # Basic damage metrics (first 3 positions)
  data <- try({
    read.table(file_path, header = TRUE, check.names = FALSE)
  }, silent = TRUE)
  
  if(inherits(data,"try-error") || nrow(data) < 3 || !"5pC>T" %in% names(data)) {
    warning(paste("Invalid data format in:", file_path))
    positions_1_3 <- c(NA, NA, NA)
  } else {
    positions_1_3 <- data[1:3, "5pC>T"]
  }
  
  # Advanced decay analysis
  decay_results <- analyze_dmg_decay(file_path)
  
  # Generate and save plot
  if (!is.na(decay_results$decay_rate)) {
    sample_name <- basename(dirname(file_path))
    plot <- try(plot_dmg_decay(decay_results, sample_name), silent = TRUE)
    
    if (!inherits(plot, "try-error") && !is.null(plot)) {
      output_dir <- dirname(dirname(file_path))
      plot_file <- file.path(output_dir, paste0(sample_name, "_dmg_decay.pdf"))
      try(ggsave(plot_file, plot, width = 8, height = 6), silent = TRUE)
    }
  }
  
  # Return combined metrics
  return(list(
    positions_1_3 = positions_1_3,
    decay_metrics = list(
      decay_rate = decay_results$decay_rate,
      N_value = decay_results$N_value,
      t_value = decay_results$t_value,
      p_value = decay_results$p_value,
      r_squared = decay_results$r_squared
    )
  ))
}

#------------------------------------------------------------------------------#
#                             Main Processing Function                         #
#------------------------------------------------------------------------------#
process_samples <- function(output_file, metadata_file = NULL) {
  # Read metadata if provided
  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
  } else {
    metadata <- NULL
  }
  
  # Get all flagstat files
  flagstat_files <- list.files(path = mapping_dir,
                               pattern = "\\.flagstat\\.log$",
                               full.names = TRUE)
  
  # Initialize results list
  results_list <- list()
  
  # Process each sample
  for (flagstat_file in flagstat_files) {
    sample_name <- sub("\\.flagstat\\.log$", "", basename(flagstat_file))
    message("Processing sample: ", sample_name)
    
    # Process flagstat file
    flagstat_results <- process_flagstat(flagstat_file)
    
    # Get duplication rate
    dup_file <- file.path(mapping_dir, paste0(sample_name, ".mapped.sorted.log"))
    dup_rate <- get_duplication_rate(dup_file)
    
    # Process MapDamage files
    sample_mapdamage_dir <- file.path(mapdamage_dir, sample_name)
    
    # Get size metrics and lambda
    lgdist_file <- file.path(sample_mapdamage_dir, "lgdistribution.txt")
    size_results <- process_lgdistribution(lgdist_file)
    
    # Get damage metrics
    ctot_file <- file.path(sample_mapdamage_dir, "5pCtoT_freq.txt")
    dmg_metrics <- extract_damage_metrics(ctot_file)
    
    # Create results row with only the requested columns
    results_list[[sample_name]] <- data.frame(
      Sample = sample_name,
      Endogenous_fraction = flagstat_results[1],
      NO_MERGED_READS = flagstat_results[2],
      DUP_rate = dup_rate,
      MEDIAN_SIZE = size_results$median_size,
      Log_Mean = size_results$log_stats[1],
      Lambda = size_results$lambda_stats[1],
      Lambda_R_squared = size_results$lambda_stats[2],
      Lambda_p_value = size_results$lambda_stats[3],
      `5P_DMG_POS1` = dmg_metrics$positions_1_3[1],
      `5P_DMG_R_squared` = dmg_metrics$decay_metrics$r_squared,
      `5P_DMG_p_value` = dmg_metrics$decay_metrics$p_value
    )
    
    # Save raw distribution data separately if needed
    if (!is.null(size_results$raw_distribution)) {
      write.table(size_results$raw_distribution,
                  file = file.path(dirname(output_file), 
                                   paste0(sample_name, "_length_dist.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  # Combine all results
  results <- bind_rows(results_list)
  
  # Merge with metadata if available
  if (!is.null(metadata)) {
    results <- left_join(results, metadata, by = "Sample")
  }
  
  # Write results
  write.table(results, 
              file = output_file, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  message("Analysis complete. Results saved in ", output_file)
  
  # Cleanup length distribution files
  cleanup_length_files(dirname(output_file))
  
  return(invisible(results))
}

#------------------------------------------------------------------------------#
#                               Cleanup Functions                              #
#------------------------------------------------------------------------------#

cleanup_length_files <- function(output_dir) {
  length_files <- list.files(path = output_dir,
                             pattern = "_length_dist.txt$", full.names = TRUE)
  if (length(length_files) > 0) {
    file.remove(length_files)
    message("Removed temporary length distribution files.")
  } else {
    message("No length distribution files found for removal.")
  }
}

cleanup_plot_files <- function(output_dir, keep_plots = TRUE) {
  if (!keep_plots) {
    plot_files <- list.files(path = output_dir,
                             pattern = "_lambda_fit.pdf$", full.names = TRUE)
    if (length(plot_files) > 0) {
      file.remove(plot_files)
      message("Removed lambda fit plot files.")
    }
  }
}

cleanup_dmg_plot_files <- function(output_dir, keep_plots = TRUE) {
  if (!keep_plots) {
    plot_files <- list.files(path = output_dir,
                             pattern = "_dmg_decay.pdf$", full.names = TRUE)
    if (length(plot_files) > 0) {
      file.remove(plot_files)
      message("Removed damage decay plot files.")
    }
  }
}

#------------------------------------------------------------------------------#
#                               Script Execution                               #
#------------------------------------------------------------------------------#

# Define output files
output_file <- "aDNA_damage_screening_MAIN.txt"
metadata_file <- "sample_metadata.txt"  # Optional, set to NULL if not available

# Run analysis
results <- process_samples(output_file, metadata_file)

# Cleanup option (set keep_plots to FALSE to remove plots)
cleanup_plot_files(dirname(output_file), keep_plots = TRUE)
cleanup_dmg_plot_files(dirname(output_file), keep_plots = TRUE)

print("Data preparation complete!")

#------------------------------------------------------------------------------#
#                 SUPPLEMENTARY MISINCORPORATIONS ANALYSIS                     #
#                   - Extract all misincorporation types                       #
#              - Calculate "other" non-deamination substitutions               #
#                  - Merge with main screening results                         #
#------------------------------------------------------------------------------#
message("\n=== Extracting misincorporations from MapDamage ===")

# =================== Function to extract all misincorporations ===============#
extract_other_misincorp <- function(file_path) {
  # Read the file
  data <- read.table(file_path, header = TRUE, comment.char = "#", 
                     stringsAsFactors = FALSE)
  
  # Filter for position 1 at both 5' and 3' ends
  pos1_5p <- data %>% filter(End == "5p", Pos == 1)
  pos1_3p <- data %>% filter(End == "3p", Pos == 1)
  
  # Calculate for 5' end (all non-deamination misincorporations)
  if (nrow(pos1_5p) > 0) {
    totals_5p <- pos1_5p %>%
      summarise(
        Total_A = sum(A),
        Total_T = sum(T),
        Total_G = sum(G),
        Total_C = sum(C),
        # Non-deamination misincorporations at 5' end
        # NOTE: R converts > to . in column names when reading tables
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
  
  # Calculate for 3' end (all non-deamination misincorporations)
  if (nrow(pos1_3p) > 0) {
    totals_3p <- pos1_3p %>%
      summarise(
        # Non-deamination misincorporations at 3' end
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

# =================== Process all samples =====================================#
all_misinc_data <- data.frame()

# Get all sample directories
sample_dirs <- list.dirs(mapdamage_dir, full.names = TRUE, recursive = FALSE)

for (sample_dir in sample_dirs) {
  sample_id <- basename(sample_dir)
  
  # File paths
  misinc_file <- file.path(sample_dir, "misincorporation.txt")
  ct_file <- file.path(sample_dir, "5pCtoT_freq.txt")
  ga_file <- file.path(sample_dir, "3pGtoA_freq.txt")
  
  # Check all files exist
  if (!file.exists(misinc_file) | !file.exists(ct_file) | !file.exists(ga_file)) {
    warning(paste("Skipping", sample_id, "- missing files"))
    next
  }
  
  # Extract 5' C>T from dedicated file (position 1)
  ct_data <- read.table(ct_file, header = TRUE)
  ct_freq <- ct_data[ct_data$pos == 1, 2]  # Second column
  
  # Extract 3' G>A from dedicated file (position 1)
  ga_data <- read.table(ga_file, header = TRUE)
  ga_freq <- ga_data[ga_data$pos == 1, 2]  # Second column
  
  # Extract all other misincorporations
  other_misinc <- extract_other_misincorp(misinc_file)
  
  if (!is.null(other_misinc)) {
    # Combine everything
    sample_data <- data.frame(
      Sample = sample_id,
      # Deamination from pre-calculated files
      `5P_DMG_POS1` = ct_freq,
      `3P_DMG_POS1` = ga_freq,
      check.names = FALSE
    )
    # Add calculated values from misincorporation.txt
    sample_data <- cbind(sample_data, other_misinc)
    
    all_misinc_data <- rbind(all_misinc_data, sample_data)
  }
}

message("Processed ", nrow(all_misinc_data), " samples for misincorporations")

# =================== Calculate "other" substitutions averages ================#
all_misinc_data <- all_misinc_data %>%
  mutate(
    other_5p_freq = (A_to_G_5p_freq + A_to_C_5p_freq + A_to_T_5p_freq +
                       T_to_C_5p_freq + T_to_G_5p_freq + T_to_A_5p_freq +
                       G_to_C_5p_freq + G_to_T_5p_freq +
                       C_to_G_5p_freq + C_to_A_5p_freq) / 10,
    
    other_3p_freq = (A_to_G_3p_freq + A_to_C_3p_freq + A_to_T_3p_freq +
                       T_to_C_3p_freq + T_to_G_3p_freq + T_to_A_3p_freq +
                       G_to_C_3p_freq + G_to_T_3p_freq +
                       C_to_G_3p_freq + C_to_A_3p_freq) / 10
  )

# Save all misincorporations
write.table(all_misinc_data, "aDNA_damage_misincorporation_and_substitutions.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)
message("Saved: aDNA_damage_misincorporation_and_substitutions.txt")

#------------------------------------------------------------------------------#
#                    STEP 4: Merge with Main Metadata                         #
#------------------------------------------------------------------------------#

message("\n=== Merging misincorporation data with main metadata ===")

# Read main metadata
metadata <- read.delim("aDNA_damage_screening_MAIN.txt")

# Check for matching IDs
misinc_ids <- unique(all_misinc_data$Sample)
metadata_ids <- unique(metadata$Sample)
matches <- sum(misinc_ids %in% metadata_ids)

message("Samples in both datasets: ", matches, " out of ", length(metadata_ids))

# Select only the columns we want to add
misinc_subset <- all_misinc_data %>%
  dplyr::select(Sample, 
                `3P_DMG_POS1`,
                other_5p_freq, 
                other_3p_freq)

# Rename columns
misinc_subset <- misinc_subset %>%
  dplyr::rename(
    `5P_other_freq` = other_5p_freq,
    `3P_other_freq` = other_3p_freq
  )

# Merge with metadata
merged_data <- metadata %>%
  left_join(misinc_subset, by = "Sample")

# Save merged data
write.table(merged_data, "aDNA_damage_screening_MAIN.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

rows_with_misinc <- sum(!is.na(merged_data$`3P_DMG_POS1`))
message("Merge complete!")
message("Rows with misincorporation data: ", rows_with_misinc, " out of ", nrow(merged_data))
message("Saved: aDNA_damage_screening_MAIN.txt")

print("Data preparation complete!")