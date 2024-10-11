

# Import packages ---------------------------------------------------------
library(plyr)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(multcomp)
library(multcompView)
library(car)

# Import files arrange data by naming -----------------------------------------
## function to import organoid data from csv files extracted by Imaris
read_plus <- function(flnm) {
  read_csv(flnm, skip = 3) %>% 
    mutate(filename = flnm)
}
## directory where the files are located
working_directory <- "/Users/majbrittbuchholz/surfdrive/Shared/PhD/Thesis/Chapter_5_Breast_on_a_Chip/Analysis/Cell coverage/Cell coverage/Data"  ## Insert here the direction to the example dataset directory
setwd(working_directory)

# import coverage
pat = ".csv"
files <- list.files(path = working_directory, pattern = pat,  recursive = TRUE, ignore.case = TRUE)
Cell_Coverage_orig <- ldply(files, read_plus)
# Create the Cell_Coverage dataframe by filtering out rows with "rotation" or "no rotation" in the Filename column
Cell_Coverage <- Cell_Coverage_orig %>%
  filter(!grepl("rotation", Filename) & !grepl("no rotation", Filename))

Cell_Coverage <- subset(Cell_Coverage, select = -c(3))
Cell_Coverage$Celltype <- Cell_Coverage$Filename
Cell_Coverage$Coating <- Cell_Coverage$Filename
Cell_Coverage$Density <- Cell_Coverage$Filename
####
Cell_Coverage$Celltype <- ifelse(grepl("HU", Cell_Coverage$Celltype), "HUVECs", Cell_Coverage$Celltype)
valid_celltypes <- c("HUVECs")
Cell_Coverage$Celltype <- ifelse(Cell_Coverage$Celltype %in% valid_celltypes, Cell_Coverage$Celltype, "MCF10A")
####
# Update the Coating column only if both "Col" and "PLL" are present in the string
Cell_Coverage$Coating <- ifelse(grepl("Col", Cell_Coverage$Coating) & grepl("PLL", Cell_Coverage$Coating), "Col1_PLL", Cell_Coverage$Coating)

# Handle the cases where only "PLL" or only "Col" are present
Cell_Coverage$Coating <- ifelse(grepl("PLL", Cell_Coverage$Coating) & !grepl("Col", Cell_Coverage$Coating), "PLL", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("Col", Cell_Coverage$Coating) & !grepl("PLL", Cell_Coverage$Coating), "Col1", Cell_Coverage$Coating)

# Handle the other specific conditions
Cell_Coverage$Coating <- ifelse(grepl("none", Cell_Coverage$Coating, ignore.case = TRUE), "None", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("LDOPA", Cell_Coverage$Coating), "LDOPA", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("Fibr", Cell_Coverage$Coating), "Fibrin", Cell_Coverage$Coating)

# Ensure that any value not in the valid coatings is set to "None"
valid_coatings <- c("PLL", "None", "Col1", "LDOPA", "Col1_PLL", "Fibrin")
Cell_Coverage$Coating <- ifelse(Cell_Coverage$Coating %in% valid_coatings, Cell_Coverage$Coating, "None")
Cell_Coverage$Density <- ifelse(grepl("7k5", Cell_Coverage$Density), "7.5", Cell_Coverage$Density)
Cell_Coverage$Density <- ifelse(grepl("15", Cell_Coverage$Density), "15", Cell_Coverage$Density)
Cell_Coverage$Density <- ifelse(grepl("30", Cell_Coverage$Density), "30", Cell_Coverage$Density)
Cell_Coverage$Density <- ifelse(grepl("60", Cell_Coverage$Density), "60", Cell_Coverage$Density)
valid_densities <- c("7.5", "15", "30", "60")

Cell_Coverage$Density <- ifelse(Cell_Coverage$Density %in% valid_densities, Cell_Coverage$Density, "30")
colnames(Cell_Coverage) <- c("Filename", "Coverage", "Celltype", "Coating", "Density")
Cell_Coverage <- Cell_Coverage %>%
  mutate(Coverage = gsub("\\;", "", Coverage)
  )
as.numeric(Cell_Coverage$Coverage)

# Coatings  ---------------------------------------------------------------
######Coatings
# Filter the data
Cell_Coverage_Coatings <- Cell_Coverage %>%
  filter(Density == "30")
Cell_Coverage_Coatings <- Cell_Coverage_Coatings %>%
  mutate(Coverage = as.numeric(Coverage)) %>%
  filter(!is.na(Coverage))

# MCF10A Coatings + Statistics --------------------------------------------
#Statistics MCF10A Coatings
# Load necessary libraries
library(dplyr)
library(tidyr)
library(multcompView)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(multcompView)

# Subset data for MCF10A cell type
Cell_Coverage_Coatings_MCF10A <- subset(Cell_Coverage_Coatings, Celltype == "MCF10A")

# Calculate summary statistics
summary_MCF10A <- Cell_Coverage_Coatings_MCF10A %>%
  group_by(Coating) %>%
  summarize(
    Average_Coverage = mean(Coverage),
    sd_Coverage = sd(Coverage),
    n = n()
  )

# Perform ANOVA
anova_result_MCF10A_Coatings <- aov(Coverage ~ Coating, data = Cell_Coverage_Coatings_MCF10A)
summary(anova_result_MCF10A_Coatings)

# Perform Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_result_MCF10A_Coatings)
# Extract p-values from Tukey's HSD test and format the data frame
p_values <- as.data.frame(tukey_result$Coating)
p_values$Comparison <- rownames(p_values)
p_values <- p_values %>%
  mutate(
    Coating1 = sub("-.*", "", Comparison),
    Coating2 = sub(".*-", "", Comparison),
    significance = case_when(
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01 ~ "**",
      `p adj` < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# Print the p-values and significance
print(p_values)

# Create a summary of comparisons and add significance
comparisons <- p_values %>%
  dplyr::select(Coating1, Coating2, `p adj`, significance)

# Prepare the data for significance labels at the bottom
significance_labels <- p_values %>%
  dplyr::select(Coating1, Coating2, significance) %>%
  mutate(ring = as.numeric(factor(Coating1, levels = c("None", "Col1", "PLL"))))

# Print the significance_labels to verify its structure
print(significance_labels)

# Merge the p-values and significance into the summary dataframe
MCF10A_Coatings <- summary_MCF10A %>%
  left_join(comparisons %>% rename(Coating = Coating1), by = "Coating") %>%
  mutate(ring = as.numeric(factor(Coating, levels = c("None", "Col1", "PLL"))))

# Print the merged summary dataframe
print(MCF10A_Coatings)

# Prepare the data for the ring chart
MCF10A_Coatings <- MCF10A_Coatings %>%
  arrange(Average_Coverage) %>%
  mutate(
    ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))),
    Average_Coverage = Average_Coverage / 100
  )


####
# Define a custom color palette with muted tones
custom_colors <- c(
  "None" = "#B2DF8A",   # Muted green
  "Col1" = "#CAB2D6",   # Muted magenta
  "PLL" = "#FB9A99",    # Muted red
  "Col1_PLL" = "#80B1D3",  # Muted cyan
  "LDOPA" = "#BC80BD",  # Muted purple
  "Fibrin" = "#8DA0CB"  # Muted blue
)

# Create a separate data frame for vertical lines and significance labels
vertical_lines <- data.frame(
  ring = c(2, 1, 3),
  label = c("NS", "*", "**"),
  y = c(0, 0, 0)
)

# Print the data frames to debug the issue
print(MCF10A_Coatings)
print(vertical_lines)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define custom colors
custom_colors <- c("None" = "#B2DF8A", "Col1" = "#CAB2D6", "PLL" = "#FB9A99")

# Basic Plot Code
MCF10A_Coatings_Graph <- ggplot(MCF10A_Coatings, aes(x = factor(ring), y = Average_Coverage, fill = Coating)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 1, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  coord_polar(theta = "y") +
  ylim(0, 1) +  # Set y-axis limits to ensure values are within the [0, 1] range
  theme_void() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Coverage by Coating for MCF10A",
    fill = "Coating"
  )

# Print the basic graph
print(MCF10A_Coatings_Graph)

# HUVECs coatings ---------------------------------------------------------
# Arrange the dataframe by increasing Average_Coverage
Cell_Coverage_Coatings_HUVECs <- subset(Cell_Coverage_Coatings, Celltype == "HUVECs")

# Calculate summary statistics
summary_HUVECs <- Cell_Coverage_Coatings_HUVECs %>%
  group_by(Coating) %>%
  summarize(
    Average_Coverage = mean(Coverage),
    sd_Coverage = sd(Coverage),
    n = n()
  )

HUVECs_Coatings <- summary_HUVECs %>%
  arrange(Average_Coverage) %>%
  mutate(
    ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))),
    Average_Coverage = Average_Coverage / 100
  )
####
# Define a custom color palette with muted tones
custom_colors <- c(
  "None" = "#B2DF8A",   # Muted green
  "Col1" = "#CAB2D6",   # Muted magenta
  "PLL" = "#FB9A99",    # Muted red
  "Col1_PLL" = "#80B1D3",  # Muted cyan
  "LDOPA" = "#BC80BD",  # Muted purple
  "Fibrin" = "#8DA0CB"  # Muted blue
)
HUVECs_Coatings <- subset(HUVECs_Coatings, Coating != "LDOPA")

# Create the stacked ring graph with statistics
HUVECs_Coatings_Graph <- ggplot(HUVECs_Coatings, aes(x = factor(ring), y = Average_Coverage, fill = Coating)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  # Add significance labels
  #geom_text(aes(label = ifelse(significance == "NS", "NS", ifelse(is.na(significance), "NA", "*")), y = 0.995), color = "black", size = 4) +
  coord_polar(theta = "y") +
  scale_y_continuous(labels = scales::percent_format()) +  # Display y-axis labels as percentages
  ylim(0, 1) +  # Set y-axis limits to ensure values are within the [0, 1] range
  theme_void() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Coverage by Coating for HUVECs",
    fill = "Coating"
  )
# Print the graph
print(HUVECs_Coatings_Graph)
#

# No stats ----------------------------------------------------------------
###with no stats
coating_order <- c("None", "Col1", "PLL", "Col1_PLL", "LDOPA", "Fibrin")
avg_Coverage_Coatings <- Cell_Coverage_Coatings %>%
  group_by(Coating, Celltype) %>%
  summarize(
    Average_Coverage = mean(Coverage),
    .groups = 'drop'  # This avoids a grouped tibble output
  )
# Reorder avg_Coverage_Density by custom order of Density
avg_Coverage_Coatings <- avg_Coverage_Coatings %>%
  mutate(Coating = factor(Coating, levels = coating_order)) %>%
  arrange(Coating)

# Prepare data
avg_Coverage_Coatings <- avg_Coverage_Coatings %>%
  group_by(Celltype) %>%
  mutate(
    ring = as.numeric(factor(Coating)))
# Convert Coating to factor to ensure proper stacking order
avg_Coverage_Coatings$Coating <- factor(avg_Coverage_Coatings$Coating)
# Convert percentages to proportions
avg_Coverage_Coatings$Average_Coverage <- avg_Coverage_Coatings$Average_Coverage / 100

########### Trying to adapt the grid lines to the length of the data
library(ggplot2)
MCF10A_Coatings <- subset(avg_Coverage_Coatings, Celltype == "MCF10A")
HUVECs_Coatings <- subset(avg_Coverage_Coatings, Celltype == "HUVECs")

# MCF10A_Coatings_Graph

# Define a custom color palette with muted tones
custom_colors <- c(
  "None" = "#B2DF8A",   # Muted green
  "Col1" = "#CAB2D6",   # Muted magenta
  "PLL" = "#FB9A99",    # Muted red
  "Col1_PLL" = "#80B1D3",  # Muted cyan
  "LDOPA" = "#BC80BD",  # Muted purple
  "Fibrin" = "#8DA0CB"  # Muted blue
)
# Arrange the dataframe by increasing Average_Coverage
MCF10A_Coatings <- MCF10A_Coatings %>%
  arrange(Average_Coverage)

# Create a new ring factor based on the order of Average_Coverage
MCF10A_Coatings <- MCF10A_Coatings %>%
  mutate(ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))))
# With percentage labels
MCF10A_Coatings_Graph <- ggplot(MCF10A_Coatings, aes(x = factor(ring), y = Average_Coverage, fill = Coating)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  # Add percentage labels
  #geom_text(aes(label = paste0(round(Average_Coverage * 100), "%"), y = 0.995), color = "black", size = 4) +
  coord_polar(theta = "y") +
  scale_y_continuous(labels = scales::percent_format()) +  # Display y-axis labels as percentages
  ylim(0, 1) +  # Set y-axis limits to ensure values are within the [0, 1] range
  theme_void() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Coverage by Coating MCF10A",
    fill = "Coating"
  )
MCF10A_Coatings_Graph

# Densities ---------------------------------------------------------------
########Density
# Filter the data
Cell_Coverage_Density <- Cell_Coverage %>%
  filter(Coating == "PLL", Celltype == "MCF10A")
Cell_Coverage_Density <- Cell_Coverage_Density %>%
  mutate(Coverage = as.numeric(Coverage)) %>%
  filter(!is.na(Coverage))
# Convert Coverage to numeric
Cell_Coverage_Density <- Cell_Coverage_Density %>%
  mutate(Coverage = as.numeric(Coverage))

Cell_Coverage_Density <- Cell_Coverage_Density %>%
  mutate(Density = as.numeric(Density))

density_order <- c("7.5", "15", "30", "60")
Cell_Coverage_Density <-Cell_Coverage_Density %>%
  mutate(Density = factor(Density, levels = density_order)) %>%
  arrange(Density)
# Calculate summary statistics
summary_Density <- Cell_Coverage_Density %>%
  group_by(Density) %>%
  summarize(
    Average_Coverage = mean(Coverage),
    sd_Coverage = sd(Coverage),
    n = n()
  )
# Perform ANOVA
anova_result_Density <- aov(Coverage ~ Density, data = Cell_Coverage_Density)
summary(anova_result_Density)
# Perform Tukey's HSD post-hoc test
tukey_result_Density <- TukeyHSD(anova_result_Density)
# Extract adjusted p-values
p_values <- tukey_result_Density$`Density`[, "p adj"]
# Create significance labels based on adjusted p-values
significance_labels <- ifelse(p_values < 0.05, "*", "NS")
# Add significance labels to the summary_Density dataframe
summary_Density <- summary_Density %>%
  mutate(significance = ifelse(row_number() <= length(significance_labels), significance_labels[row_number()], NA))

# Prepare the data for the ring chart
MCF10A_Density <- summary_Density %>%
  arrange(Density) %>%
  mutate(
    ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))),
    Average_Coverage = Average_Coverage / 100
  )
####
# Define a custom color palette with muted tones
custom_colors_density <- c(
  "7.5" = "#B2DF8A",   # Muted green
  "15" = "#CAB2D6",   # Muted magenta
  "30" = "#FB9A99",    # Muted red
  "60" = "#80B1D3"  # Muted cyan
)
# Create the stacked ring graph with statistics
MCF10A_Density_Graph <- ggplot(MCF10A_Density, aes(x = factor(ring), y = Average_Coverage, fill = Density)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  scale_fill_manual(values = custom_colors_density) +  # Apply custom colors
  # Add significance labels
  coord_polar(theta = "y") +
  scale_y_continuous(labels = scales::percent_format()) +  # Display y-axis labels as percentages
  ylim(0, 1) +  # Set y-axis limits to ensure values are within the [0, 1] range
  theme_void() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Coverage by Density for MCF10A",
    fill = "Density"
  )

# Print the graph
print(MCF10A_Density_Graph)

# no statistics density ---------------------------------------------------
avg_Coverage_Density <- Cell_Coverage_Density %>%
  group_by(Density) %>%
  summarize(Average_Coverage = mean(Coverage, na.rm = TRUE)/100)


avg_Coverage_Density <- avg_Coverage_Density %>%
  mutate(
    ring = as.numeric(factor(Density, levels = density_order))
  )

avg_Coverage_Density <- avg_Coverage_Density %>%
  mutate(
    Average_Coverage = as.numeric(Average_Coverage)
  )

# Arrange the dataframe by increasing Average_Coverage
avg_Coverage_Density <- avg_Coverage_Density %>%
  arrange(Average_Coverage)

# Create a new ring factor based on the order of Average_Coverage
avg_Coverage_Density <- avg_Coverage_Density %>%
  mutate(ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))))

# # Create the stacked ring chart
# ggplot(avg_Coverage_Density, aes(x = factor(ring), y = Average_Coverage, fill = reorder(factor(Density), density_order))) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   scale_y_continuous(labels = scales::percent_format()) +  # Display y-axis labels as percentages
#   ylim(0, 1) +  # Set y-axis limits to ensure values are within the [0, 1] range
#   theme_void() +
#   theme(legend.position = "bottom") +
#   labs(
#     title = "Average Coverage by Density for Coating = PLL and Celltype = MCF10A",
#     fill = "Density"
#   )
# 
# Rotation vs no rotation -------------------------------------------------
# Create the Flipping dataframe by filtering rows with "rotation" or "no rotation" in the Filename column
Flipping <- Cell_Coverage_orig %>%
  filter(grepl("rotation", Filename) | grepl("norotation", Filename)| grepl("No_Flip", Filename)| grepl("No_Coat", Filename))
#
Flipping$Filename <- ifelse(grepl("norotation", Flipping$Filename), "No Rotation", Flipping$Filename)
Flipping$Filename <- ifelse(grepl("No_Flip", Flipping$Filename), "No Rotation", Flipping$Filename)
valid <- c("No Rotation")

Flipping$Filename <- ifelse(Flipping$Filename %in% valid, Flipping$Filename, "Rotation")
# Arrange by Coverage and create a ring variable
Flipping <- Flipping %>%
  arrange(Coverage) %>%
  mutate(
    ring = as.numeric(factor(Coverage, levels = unique(Coverage))),
    Coverage = Coverage / 100
  )

# Assuming Flipping dataframe is already correctly filtered
# Create the custom color scheme
custom_colors_flipping <- c("No Rotation" = "#7fc97f", "Rotation" = "#beaed4")

# Generate the ring chart
Flipping_Graph <- ggplot(Flipping, aes(x = factor(Filename), y = Coverage, fill = Filename)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1.05, alpha = 0.5, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  # Add percentage labels
  #geom_text(aes(label = paste0(round(Coverage * 100), "%"), y = 0.995), color = "black", size = 4) +
  # Apply custom colors
  scale_fill_manual(values = custom_colors_flipping) +
  # Polar coordinates for the ring effect
  coord_polar(theta = "y") +
  # y-axis labels as percentages
  scale_y_continuous(labels = scales::percent_format()) +
  # Ensure values are within the [0, 1] range
  ylim(0, 1) +
  # Simplified theme for clarity
  theme_void() +
  # Customize legend position
  theme(legend.position = "bottom") +
  # Title and axis labels
  labs(
    title = "Coverage by Filename",
    fill = "Filename"
  )

# Display the graph
Flipping_Graph

###
# Arrange the dataframe by increasing Average_Coverage
Cell_Coverage_None <- subset(Cell_Coverage, Coating == "None")
# Convert Coverage column to numeric
Cell_Coverage_None$Coverage <- as.numeric(Cell_Coverage_None$Coverage)
summary_None <- Cell_Coverage_None %>%
  summarise(
    Average_Coverage = mean(Coverage),
    sd_Coverage = sd(Coverage),
    n = n()
  )




