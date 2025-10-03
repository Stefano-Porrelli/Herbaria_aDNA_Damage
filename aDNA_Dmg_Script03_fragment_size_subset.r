#!/usr/bin/env Rscript
#------------------------------------------------------------------------------#
#            Patterns of aDNA Damage Through Time end Environments             #
#                      – lessons from herbarium specimens  -                   #
#                                                                              #
#                                   Script 03                                  #
#                                                                              #
#                      DATA ANALYSIS III - Fragment Size Subset:               #
#                                                                              #
# Uses a subset of the data, which contains information of peak sizes of gDNA  #
#            generated with capillary electrophoresis (TapeStation)            #
#------------------------------------------------------------------------------#

# Install required packages
# install.packages(c("gridExtra", "grid", "ggplot2"))

# Load required libraries
library(ggplot2)
library(grid)
library(gridExtra)

# Read data
d <- read.delim("aDNA_damage_screening_SUBSET_TAPESTATION_gDNA.txt", 
                row.names = NULL, stringsAsFactors = TRUE)

# Clean the data - convert "-" to NA and make numeric
d$gDNA_peak_size <- as.numeric(as.character(d$gDNA_peak_size))
d$Library_Peak_size <- as.numeric(as.character(d$Library_Peak_size))
d$MEDIAN_SIZE <- as.numeric(as.character(d$MEDIAN_SIZE))

# Define color palette
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")

#------------------------------------------------------------------------------#
#     1 - Regression analysis: Tapestation peaks vs merged fragment sizes      #
#------------------------------------------------------------------------------#

# Calculate correlations
# Plot 1: gDNA all samples
model1 <- lm(gDNA_peak_size ~ MEDIAN_SIZE, data = d)
r2_1 <- summary(model1)$r.squared
p_1 <- summary(model1)$coefficients[2, 4]

# Plot 2: Library all samples  
model2 <- lm(Library_Peak_size ~ MEDIAN_SIZE, data = d)
r2_2 <- summary(model2)$r.squared
p_2 <- summary(model2)$coefficients[2, 4]

# Plot 3: gDNA peak size vs Library peak size
model3 <- lm(Library_Peak_size ~ gDNA_peak_size, data = d)
r2_3 <- summary(model3)$r.squared
p_3 <- summary(model3)$coefficients[2, 4]

# Plot 1
P1 <- ggplot(d, aes(x = MEDIAN_SIZE, y = gDNA_peak_size, color = Genus)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Median Fragment Size - Merged (bp)",
         y = "gDNA Peak Size (bp)",
         color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    annotate("text", size = 3.5, x = max(d$MEDIAN_SIZE, na.rm = TRUE) - 2, 
             y = min(d$gDNA_peak_size, na.rm = TRUE) + 40,
             label = sprintf("italic(R)^2 == %.3f~~italic(p) == %.3g", r2_1, p_1),
             parse = TRUE, hjust = 1)

# Plot 2
P2 <- ggplot(d, aes(x = MEDIAN_SIZE, y = Library_Peak_size, color = Genus)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Median Fragment Size - Merged (bp)",
         y = "Library Peak Size (bp)",
         color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    annotate("text", size = 3.5, x = max(d$MEDIAN_SIZE, na.rm = TRUE) - 2, 
             y = 175,
             label = sprintf("italic(R)^2 == %.3f~~italic(p) == %.3g", r2_2, p_2),
             parse = TRUE, hjust = 1)

# Plot 3: gDNA peak size vs Library peak size
P3 <- ggplot(d, aes(x = gDNA_peak_size, y = Library_Peak_size, color = Genus)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "gDNA Peak Size (bp)",
         y = "Library Peak Size (bp)",
         color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    annotate("text", size = 3.5, x = 300, 
             y =175,
             label = sprintf("italic(R)^2 == %.3f~~italic(p) == %.3g", r2_3, p_3),
             parse = TRUE, hjust = 1)

# Remove legends from plots and add labels
P1 <- P1 + theme(legend.position = "none") + labs(tag = "a")
P2 <- P2 + theme(legend.position = "none") + labs(tag = "b")
P3 <- P3 + theme(legend.position = "none") + labs(tag = "c")

# Extract legend from one plot
legend <- ggplot(d, aes(x = MEDIAN_SIZE, y = Library_Peak_size, color = Genus)) +
    geom_point() +
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = NULL))

legend_grob <- ggplotGrob(legend)$grobs[[which(sapply(ggplotGrob(legend)$grobs, function(x) x$name) == "guide-box")]]

# Combine plots: P1 and P2 on top, P3 centered at bottom (same size)
top_row <- arrangeGrob(P1, P2, ncol = 2)
bottom_row <- arrangeGrob(nullGrob(), P3, nullGrob(), ncol = 3, widths = c(1, 2, 1))

# Combine all plots with legend
grid.arrange(top_row, bottom_row, legend_grob, 
             heights = c(5, 5, 1), nrow = 3)


#------------------------------------------------------------------------------#
#          2 - Regression analysis: gDNA peak size ~ collection year           #
#------------------------------------------------------------------------------#

# gDNA peak size vs Collection_year - Combined regression
model_year_combined <- lm(gDNA_peak_size ~ Collection_year, data = d)
r2_year_combined <- summary(model_year_combined)$r.squared
p_year_combined <- summary(model_year_combined)$coefficients[2, 4]

# Plot
ggplot(d, aes(x = Collection_year, y = gDNA_peak_size, color = Genus)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +  # Single regression line in black
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Collection Year",
         y = "gDNA Peak Size (bp)",
         color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "bottom") +
    annotate("text", size = 3.5, x = max(d$Collection_year, na.rm = TRUE) - 10,
             y = 70,
             label = sprintf("italic(R)^2 == %.3f~~italic(p) == %.3g", r2_year_combined, p_year_combined),
             parse = TRUE, hjust = 1)

# Exclude outliers
d_filtered <- d[!d$Sample %in% c("HV0061", "HV0081"), ]

# gDNA peak size vs Collection_year - Combined regression (excluding outliers)
model_year_combined_filtered <- lm(gDNA_peak_size ~ Collection_year, data = d_filtered)
r2_year_combined_filtered <- summary(model_year_combined_filtered)$r.squared
p_year_combined_filtered <- summary(model_year_combined_filtered)$coefficients[2, 4]

# Plot
ggplot(d_filtered, aes(x = Collection_year, y = gDNA_peak_size, color = Genus)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +  # Single regression line in black
    scale_color_manual(values = genus_colors,
                       labels = c(expression(italic("Hordeum")),
                                  expression(italic("Oryza")))) +
    labs(x = "Collection Year",
         y = "gDNA Peak Size (bp)",
         color = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), legend.position = "bottom") +
    annotate("text", size = 3.5, x = max(d_filtered$Collection_year, na.rm = TRUE) - 10,
             y = 50,
             label = sprintf("italic(R)^2 == %.3f~~italic(p) == %.3g", r2_year_combined_filtered, p_year_combined_filtered),
             parse = TRUE, hjust = 1)
