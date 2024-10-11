

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
## function to import csv file
read_plus <- function(flnm) {
  read_csv(flnm, skip = 3) %>% #For our data, column headers start in row 4; adjust if needed
    mutate(filename = flnm)
}
## directory where the summary csv file is located. Fill in your file path
working_directory <- "/YourFilePath"  ## Insert here the direction to the example dataset directory
setwd(working_directory)

# import coverage data
pat = ".csv"
files <- list.files(path = working_directory, pattern = pat,  recursive = TRUE, ignore.case = TRUE)
Cell_Coverage_orig <- ldply(files, read_plus)
# Create the Cell_Coverage dataframe by filtering out rows with "rotation" or "no rotation" in the Filename column. we only want to look at data for rotation
Cell_Coverage <- Cell_Coverage_orig %>%
  filter(!grepl("rotation", Filename) & !grepl("no rotation", Filename))

#Remove column containing the summary file name
Cell_Coverage <- subset(Cell_Coverage, select = -c(3))

#Add three columns containing the Filename
Cell_Coverage$Celltype <- Cell_Coverage$Filename
Cell_Coverage$Coating <- Cell_Coverage$Filename
Cell_Coverage$Density <- Cell_Coverage$Filename

#### Identify celltype based on pattern in Filename
Cell_Coverage$Celltype <- ifelse(grepl("HU", Cell_Coverage$Celltype), "HUVECs", Cell_Coverage$Celltype)
valid_celltypes <- c("HUVECs")
Cell_Coverage$Celltype <- ifelse(Cell_Coverage$Celltype %in% valid_celltypes, Cell_Coverage$Celltype, "MCF10A")

#### Identify the coating based on patterns in the coating column
# Update the Coating column only if both "Col" and "PLL" are present in the string
Cell_Coverage$Coating <- ifelse(grepl("Col", Cell_Coverage$Coating) & grepl("PLL", Cell_Coverage$Coating), "Col1_PLL", Cell_Coverage$Coating)

# Handle the cases where only "PLL" or only "Col" are present
Cell_Coverage$Coating <- ifelse(grepl("PLL", Cell_Coverage$Coating) & !grepl("Col", Cell_Coverage$Coating), "PLL", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("Col", Cell_Coverage$Coating) & !grepl("PLL", Cell_Coverage$Coating), "Col1", Cell_Coverage$Coating)

# Handle the other specific conditions
Cell_Coverage$Coating <- ifelse(grepl("none|None", Cell_Coverage$Coating, ignore.case = TRUE), "None", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("LDOPA", Cell_Coverage$Coating), "LDOPA", Cell_Coverage$Coating)
Cell_Coverage$Coating <- ifelse(grepl("Fibr", Cell_Coverage$Coating), "Fibrin", Cell_Coverage$Coating)

# Ensure that any value not in the valid coatings is set to "None"
valid_coatings <- c("PLL", "None", "Col1", "LDOPA", "Col1_PLL", "Fibrin")
Cell_Coverage$Coating <- ifelse(Cell_Coverage$Coating %in% valid_coatings, Cell_Coverage$Coating, "None")

# Identify seeding density based on pattern in Filenmae
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
# Filter the data by same seeding density but different coatings; adjust if needed
Cell_Coverage_Coatings <- Cell_Coverage %>%
  filter(Density == "30")
Cell_Coverage_Coatings <- Cell_Coverage_Coatings %>%
  mutate(Coverage = as.numeric(Coverage)) %>%
  filter(!is.na(Coverage))

# MCF10A Coatings + Statistics --------------------------------------------
#Statistics MCF10A Coatings
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
p_values_MCF10A_Coatings <- as.data.frame(tukey_result$Coating)
p_values_MCF10A_Coatings$Comparison <- rownames(p_values_MCF10A_Coatings)
p_values_MCF10A_Coatings <- p_values_MCF10A_Coatings %>%
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
print(p_values_MCF10A_Coatings)

# Create a summary of comparisons and add significance
comparisons_MCF10A <- p_values_MCF10A_Coatings %>%
  dplyr::select(Coating1, Coating2, `p adj`, significance)

# Prepare the data for significance labels at the bottom
significance_labels <- p_values_MCF10A_Coatings %>%
  dplyr::select(Coating1, Coating2, significance) %>%
  mutate(ring = as.numeric(factor(Coating1, levels = c("None", "Col1", "PLL"))))

# Print the significance_labels to verify its structure
print(significance_labels)

# Merge the p-values and significance into the summary dataframe
MCF10A_Coatings <- summary_MCF10A %>%
  left_join(comparisons_MCF10A %>% rename(Coating = Coating1), by = "Coating") %>%
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

# HUVECs coatings + Statistics ---------------------------------------------------------
#Statistics HUVEC Coatings
# Subset data for MCF10A cell type
Cell_Coverage_Coatings_HUVECs <- subset(Cell_Coverage_Coatings, Celltype == "HUVECs")

# Calculate summary statistics
summary_HUVECs <- Cell_Coverage_Coatings_HUVECs %>%
  group_by(Coating) %>%
  summarize(
    Average_Coverage = mean(Coverage),
    sd_Coverage = sd(Coverage),
    n = n()
  )

# Perform ANOVA
anova_result_HUVECs_Coatings <- aov(Coverage ~ Coating, data = Cell_Coverage_Coatings_HUVECs)
summary(anova_result_HUVECs_Coatings)

# Perform Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_result_HUVECs_Coatings)
# Extract p-values from Tukey's HSD test and format the data frame
p_values_HUVECs_Coatings <- as.data.frame(tukey_result$Coating)
p_values_HUVECs_Coatings$Comparison <- rownames(p_values_HUVECs_Coatings)
p_values_HUVECs_Coatings <- p_values_HUVECs_Coatings %>%
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
print(p_values_HUVECs_Coatings)

# Create a summary of comparisons and add significance
comparisons_HUVECs <- p_values_HUVECs_Coatings %>%
  dplyr::select(Coating1, Coating2, `p adj`, significance)

# Prepare the data for significance labels at the bottom
significance_labels <- p_values_HUVECs_Coatings %>%
  dplyr::select(Coating1, Coating2, significance) %>%
  mutate(ring = as.numeric(factor(Coating1, levels = c("None", "Col1", "PLL", "Col1_PLL", "Fibrin"))))

# Print the significance_labels to verify its structure
print(significance_labels)

# Merge the p-values and significance into the summary dataframe
HUVECs_Coatings <- summary_HUVECs %>%
  left_join(comparisons_HUVECs %>% rename(Coating = Coating1), by = "Coating") %>%
  mutate(ring = as.numeric(factor(Coating, levels = c("None", "Col1", "PLL", "Col1_PLL", "Fibrin"))))

# Print the merged summary dataframe
print(HUVECs_Coatings)

# Prepare the data for the ring chart
HUVECs_Coatings <- HUVECs_Coatings %>%
  arrange(Average_Coverage) %>%
  mutate(
    ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))),
    Average_Coverage = Average_Coverage / 100
  )

# Make the graphs without statistics indications (add manually later) ---------
# Define a custom color palette with muted tones
custom_colors <- c(
  "None" = "#B2DF8A",   # Muted green
  "Col1" = "#CAB2D6",   # Muted magenta
  "PLL" = "#FB9A99",    # Muted red
  "Col1_PLL" = "#80B1D3",  # Muted cyan
  "LDOPA" = "#BC80BD",  # Muted purple
  "Fibrin" = "#8DA0CB"  # Muted blue
)

# MCF10A_Coatings_Graph
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
  geom_text(aes(label = paste0(round(Average_Coverage * 100), "%"), y = 0.995), color = "black", size = 4) +
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

# HUVECs_Coatings_Graph
# Arrange the dataframe by increasing Average_Coverage
HUVECs_Coatings <- HUVECs_Coatings %>%
  arrange(Average_Coverage)

# Create a new ring factor based on the order of Average_Coverage
HUVCECs_Coatings <- HUVECs_Coatings %>%
  mutate(ring = as.numeric(factor(Average_Coverage, levels = unique(Average_Coverage))))
# With percentage labels
HUVECs_Coatings_Graph <- ggplot(HUVECs_Coatings, aes(x = factor(ring), y = Average_Coverage, fill = Coating)) +
  # Background layer for the full ring
  geom_bar(aes(y = 1), fill = "lightgrey", color = NA, width = 1, stat = "identity") +
  # Data layer
  geom_bar(stat = "identity", width = 1) +
  # Grid lines
  geom_hline(yintercept = seq(0, 0.8, by = 0.2), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(1), color = "darkgrey", linetype = "dashed") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  # Add percentage labels
  geom_text(aes(label = paste0(round(Average_Coverage * 100), "%"), y = 0.995), color = "black", size = 4) +
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
HUVECs_Coatings_Graph
# Densities ---------------------------------------------------------------
# Filter the data; we only tried different densities on MCF10A with PLL coating. Adjust if needed
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
  # Add percentage labels
  geom_text(aes(label = paste0(round(Average_Coverage * 100), "%"), y = 0.995), color = "black", size = 4) +
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








