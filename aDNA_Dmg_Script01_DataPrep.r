#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
#                 aDNA Damage through time -  Script 01                        #
#                             Data preparation                                 #
#------------------------------------------------------------------------------#
# Install required packages
install.packages(c("dplyr", "ggplot2", "stringr", "colorspace",
                   "tidyr", "ggrepel", "ggtext", "sf", "geodata",
                   "terra", "purrr", "maps"))
# Load libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(colorspace)
library(tidyr)
library(ggrepel)
library(ggtext)
library(sf)
library(geodata)
library(terra)
library(purrr)
library(maps)

# Read data matrix in R
d <- read.delim("aDNA_damage_screening.txt",
                row.names = 1, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------#
#               1 - Descriptive statistic/data screening                       #
#------------------------------------------------------------------------------#
# Create table to count no. of accession for each species, divided by Herbarium
# Count total number of accessions for each herbarium
species_summary <- table(paste(d$Genus, d$Species),
                         d$Herbarium)
# Convert to data frame and add total count
species_summary_df <- as.data.frame.matrix(species_summary)
species_summary_df$Total <- rowSums(species_summary_df)
# Move Total column to the front
species_summary_df <- species_summary_df [, c("Total",
                                              setdiff(names(species_summary_df),
                                                      "Total"))]
# Add Species column
species_summary_df$Species <- rownames(species_summary_df)
rownames(species_summary_df) <- NULL
# Reorder columns to put species first
species_summary_df <- species_summary_df[, c("Species", "Total",
                                             setdiff(names(species_summary_df),
                                                     c("Species", "Total")))]
# Add ranges of collection year for each species divided by Herbarium
# Calculate the age range for each unique species
age_range_df <- d %>%
  group_by(Genus, Species) %>%
  summarise(Oldest_year = min(Collection_year, na.rm = TRUE),
            Youngest_year = max(Collection_year, na.rm = TRUE)) %>%
  mutate(Age_range = paste0(Oldest_year, "-", Youngest_year)) %>%
  ungroup()
# Create a unique identifier for merging
age_range_df$Species <- paste(age_range_df$Genus, age_range_df$Species)
# Merge the species summary with the age range data frame
final_species_summary <- species_summary_df %>%
  left_join(age_range_df[, c("Species", "Age_range")], by = "Species")
# Print the final results
print(final_species_summary, row.names = FALSE)

#------------------------------------------------------------------------------#
#                   define colours for species                                 #
#------------------------------------------------------------------------------#
# Get unique genera and create base colors for each
genus_colors <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
# Define the adjust_color function
adjust_color <- function(hex_color, factor) {
  # Convert hex to RGB
  rgb_col <- col2rgb(hex_color)
  # Adjust each channel
  adjusted <- rgb_col * factor
  # Ensure values stay within 0-255
  adjusted <- pmin(255, pmax(0, adjusted))
  # Convert back to hex
  return(rgb(adjusted[1] / 255, adjusted[2] / 255, adjusted[3] / 255))
}
# Create a combined Genus-Species dataframe from all possible sources
species_df <- unique(data.frame(
  Genus = sapply(strsplit(as.character(species_summary_df$Species), " "), `[`, 1),
  Species = species_summary_df$Species
))
# Initialize color map
species_colors <- list()
# Loop over each genus and assign adjusted colors per species
for (genus in unique(species_df$Genus)) {
  base_color <- genus_colors[genus]
  species_subset <- species_df[species_df$Genus == genus, ]
  n_species <- nrow(species_subset)
  factors <- seq(0.7, 1.3, length.out = n_species)  # slightly wider range
  for (i in seq_len(n_species)) {
    species_name <- species_subset$Species[i]
    species_colors[[species_name]] <- adjust_color(base_color, factors[i])
  }
}
# Convert to named vector
species_colors <- unlist(species_colors)

#------------------------------------------------------------------------------#
#           2  - Plot N. of accession/species by Herbarium                     #
#------------------------------------------------------------------------------#
# Reshape the data to have Herbarium on x-axis
plot_data <- tidyr::pivot_longer(species_summary_df,
                                 cols = -c(Species, Total),
                                 names_to = "Herbarium",
                                 values_to = "Count")
# Add genus information to plot_data
plot_data$Genus <- sapply(strsplit(as.character(plot_data$Species), " "), function(x) x[1])
# Summarize total counts per herbarium
herbarium_order <- plot_data %>%
    dplyr::group_by(Herbarium) %>%
    dplyr::summarise(total_count = sum(Count)) %>%
    dplyr::arrange(desc(total_count)) %>%
    dplyr::pull(Herbarium)
# Reorder Herbarium factor levels
plot_data$Herbarium <- factor(plot_data$Herbarium, levels = herbarium_order)

# Create the plot
P1_accessions <- ggplot(plot_data,
                                  aes(x = Herbarium,
                                      y = Count,
                                      fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(),
          legend.position = "bottom",
          legend.text = element_text(face = "italic")) +
    labs(title = NULL,
         x = "Herbarium",
         y = "No. of Herbarium samples",
         fill = NULL) +
    scale_fill_manual(values = species_colors)
# print plot
P1_accessions

#------------------------------------------------------------------------------#
#           3 - Histogram of collection year by genera                         #
#------------------------------------------------------------------------------#
d$Species <- paste(d$Genus, d$Species)
P2_temporal_all <- ggplot(d,
       aes(x = Collection_year,
           fill = Species)) +
    geom_histogram(binwidth = 5, position = "stack", color = "white", size = 0.2) +
    scale_fill_manual(values = species_colors) +
    labs(title = NULL,
         x = "Collection Year (5-years bins)",
         y = "No. of Herbarium samples",
         fill = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(face = "italic"))
# print plot
P2_temporal_all

#------------------------------------------------------------------------------#
#                       4 - Prepare data for Maps                              #
#------------------------------------------------------------------------------#
# Create a data frame of unique countries and their counts
country_counts <- d %>%
  group_by(Country, Genus) %>%
  summarise(Count = n(), .groups = "drop")
# Get map data
world_map <- map_data("world") %>%
  filter(lat > -60)  # Remove Antarctica
# Check for unmatched countries
unmatched_countries <- country_counts %>%
  anti_join(world_map, by = c("Country" = "region")) %>%
  pull(Country) %>%
  unique()
# Make sure all countries are included
if (length(unmatched_countries) > 0) {
  warning("The following countries are not recognized: ",
          paste(unmatched_countries, collapse = ", "))
}
# Calculate country centroids for plotting
country_centroids <- world_map %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat), .groups = "drop") %>%
  left_join(country_counts, by = c("region" = "Country")) %>%
  mutate(Count = replace_na(Count, 0))
# Filter out countries with a count of 0
country_centroids <- country_centroids %>% filter(Count > 0)
country_centroids <- country_centroids %>% mutate(Country = region)
# Define custom size breaks and labels
size_breaks <- c(10, 20, 30, 40)
# Define colors for the two genera
color_map <- c("Hordeum" = "#663399", "Oryza" = "#00A878")
# Create separate datasets for each genus
hordeum_data <- country_centroids %>% filter(Genus == "Hordeum") %>%
  mutate(Genus = "Hordeum")
oryza_data <- country_centroids %>% filter(Genus == "Oryza") %>%
  mutate(Genus = "Oryza")

#------------------------------------------------------------------------------#
#              5 - Map of Geographical Distrubution Hordeum/Oryza              #
#------------------------------------------------------------------------------#
# Plot datapoints on the map
P3_Geography <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "lightgrey", color = NA) +
  geom_point(data = hordeum_data,
             aes(x = long, y = lat, size = Count, color = "Hordeum"),
             alpha = 0.25) +
  geom_point(data = oryza_data,
             aes(x = long, y = lat, size = Count, color = "Oryza"),
             alpha = 0.25) +
  scale_size_continuous(breaks = size_breaks,
                        labels = size_breaks,
                        range = c(1, 15),
                        name = "Sample Count") +
  scale_color_manual(values = color_map, name = NULL) +
  theme_minimal() +
  coord_fixed(1.3) +
  labs(title = NULL) +
  theme(panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_markdown(),
    legend.text = element_markdown()
  ) # +
# Add countries labels (optional), uncomment if needed
#   geom_text_repel(
#   data = country_centroids,
#   aes(x = long, y = lat, label = paste(Country, "(n =", Count, ")")),
#   size = 3,
#   segment.color = "grey50",
#   segment.alpha = 0.5,
#   segment.size = 0.5,
#   max.overlaps = Inf,
#   box.padding = 0.3,
#   point.padding = 0.3,
#   force = 15,
#   max.iter = 20000,
#   direction = "both"
# )
# Plot
P3_Geography

#------------------------------------------------------------------------------#
#  6 - Map of Geographical Distrubution Hordeum/Oryza + mean annual temp       #
#------------------------------------------------------------------------------#
# Download and load WorldClim temperature data
temp <- worldclim_global(var = "tavg", res = 10, path = tempdir())
annual_temp <- mean(temp) # Calculate annual mean from monthly averages
# Convert temperature to a data frame for ggplot
temp_points <- as.data.frame(annual_temp, xy=TRUE)
colnames(temp_points) <- c("long", "lat", "temp")
# Remove Antarctica (comment/uncomment as needed)
temp_points <- temp_points %>% filter(lat > -60)

# Plot datapoints on map
P4_MapTemp <- ggplot() +
  # Add temperature layer
  geom_tile(data = temp_points,
            aes(x = long, y = lat, fill = temp)) +
  scale_fill_gradientn(colors = c("#313695",
                                  "#4575b4",
                                  "#74add1",
                                  "#abd9e9",
                                  "#e0f3f8",
                                  "#ffffbf",
                                  "#fee090",
                                  "#fdae61",
                                  "#f46d43",
                                  "#d73027",
                                  "#a50026"),
                       name = "Mean Annual\nTemperature (°C)",
                       breaks = seq(-30, 30, by = 10),
                       limits = c(-31, 31),
                       guide = guide_colorbar(direction = "horizontal",
                                              title.position = "top",
                                              barwidth = 10,
                                              barheight = 0.5)) +
  # Add points for species data - NOW WITH HOLLOW POINTS
  geom_point(data = hordeum_data,
             aes(x = long, y = lat, size = Count, color = "Hordeum"),
             fill = NA,  # This makes the points hollow
             shape = 21, # Circle with both fill and color
             stroke = 1, # Thicker outline
             alpha = 1) + # Increased alpha for better visibility
  geom_point(data = oryza_data,
             aes(x = long, y = lat, size = Count, color = "Oryza"),
             fill = NA,  # This makes the points hollow
             shape = 21, # Circle with both fill and color
             stroke = 1, # Thicker outline
             alpha = 1) + # Increased alpha for better visibility
  scale_size_continuous(breaks = size_breaks,
                        labels = size_breaks,
                        range = c(1, 15),
                        name = "Sample Count") +
  scale_color_manual(values = color_map,
                     name = NULL,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  theme_minimal() +
  coord_fixed(1.3) +
  labs(title = NULL,
       subtitle = NULL) +
  guides(fill = guide_colorbar(direction = "horizontal", order = 1,
                               title.position = "top",
                               barwidth = 10, barheight = 0.5),
    color = guide_legend(order = 2),    # Legend for genera
    size = guide_legend(order = 3)      # Legend for Sample Count
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm")  # Adjust spacing for readability
  )

# Plot
P4_MapTemp

#------------------------------------------------------------------------------#
#  7 - Map of Geographical Distrubution Hordeum/Oryza + annual precipitation   #
#------------------------------------------------------------------------------#
# Download and load WorldClim precipitation data
precip <- worldclim_global(var = "prec", res = 10, path = tempdir())
annual_precip <- sum(precip) # Calculate annual prepicitation
# Convert precipitation to a data frame for ggplot
precip_points <- as.data.frame(annual_precip, xy=TRUE)
colnames(precip_points) <- c("long", "lat", "precip")
# Remove Antarctica (comment/uncomment as needed)
precip_points <- precip_points %>% filter(lat > -60)  # Remove Antarctica

# Create the plot with precipitation data
P5_MapPrecipit <- ggplot() +
  # Add precipitation layer
  geom_tile(data = precip_points,
            aes(x = long, y = lat, fill = precip)) +
  scale_fill_gradientn(colors = c("#f7fbff",
                                  "#9ecae1",
                                  "#6baed6",
                                  "#4292c6",
                                  "#2171b5",
                                  "#08519c",
                                  "#08306b"),
                       name = "Annual Precipitation (mm)",
                       breaks = seq(0, 5000, by = 1000),
                       limits = c(0, 5000),
                       guide = guide_colorbar(direction = "horizontal",
                                              title.position = "top",
                                              barwidth = 10,
                                              barheight = 0.5)) +
  # Add points for species data - NOW WITH HOLLOW POINTS
  geom_point(data = hordeum_data,
             aes(x = long, y = lat, size = Count, color = "Hordeum"),
             fill = NA,  # This makes the points hollow
             shape = 21, # Circle with both fill and color
             stroke = 1, # Thicker outline
             alpha = 1) + # Increased alpha for better visibility
  geom_point(data = oryza_data,
             aes(x = long, y = lat, size = Count, color = "Oryza"),
             fill = NA,  # This makes the points hollow
             shape = 21, # Circle with both fill and color
             stroke = 1, # Thicker outline
             alpha = 1) + # Increased alpha for better visibility
  scale_size_continuous(breaks = size_breaks,
                        labels = size_breaks,
                        range = c(1, 15),
                        name = "Sample Count") +
  scale_color_manual(values = color_map,
                     name = NULL,
                     labels = c(expression(italic("Hordeum")),
                                expression(italic("Oryza")))) +
  theme_minimal() +
  coord_fixed(1.3) +
  labs(title = NULL,
       subtitle = NULL) +
  guides(fill = guide_colorbar(direction = "horizontal", order = 1,
                               title.position = "top", barwidth = 10,
                               barheight = 0.5),
    color = guide_legend(order = 2),    # Legend for Hordeum/Oryza
    size = guide_legend(order = 3)      # Legend for Sample Count
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm")  # Adjust spacing for readability
  )

# Plot
P5_MapPrecipit
